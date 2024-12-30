#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 19:05:36 2024

@author: jiachangfu
"""


#from collections import Counter
import pysam



def extract_reads_supporting_snp(bamfile, chrom, pos, target_base, invert=False ):
    """
    Parameters : target_base 
    
    The function is aming to extract all the reads in the given genomic region with the target genotype.
    1)Pysam pileup to extract the reads in the give position without any quality filter in this step. We want reads as much as possibile and filter them in the next steps.
    2)For this position, we extracted all the reads's genotype. The reads have the target genotype would be served as target reads.
    3)Byyond SNPs, we also labelled indel, deletion in the target reads or not based on is_del, indel value.
    
    DEBUG adding: 4) Some reads have clip mapping which refer to only partial reads can be mapped on this region, while the rest of them were not. This can
    cause some bugs, we just delete them through cigartuples.
    """
    
    target_reads = []
    non_target_reads = []
    for pileupcolumn in bamfile.pileup(chrom, pos, pos+1, truncate=True,min_base_quality=0, min_mapping_quality=0,stepper="nofilter"):
        if pileupcolumn.reference_pos != pos:
            continue
        for pileupread in pileupcolumn.pileups:
            
            #Soft clipping delete:  This is soft clipping sign of cigar [136M15S];  This is for none soft clipping [117M13I21M]. Only match and Indel.
            cigartuples = pileupread.alignment.cigartuples
            if cigartuples is not None:
                has_soft_clipping = any(op == 4 for op, length in cigartuples)
                if has_soft_clipping:
                    continue  # 跳过含有 soft clipping 的 reads

            #Separate deletion, Indel 
            if pileupread.is_del:
                if target_base == "-":
                    target_reads.append(pileupread.alignment)
                else:
                    non_target_reads.append(pileupread.alignment)
                continue
            elif pileupread.indel < 0:
                if target_base == "-":
                    target_reads.append(pileupread.alignment)
                else:
                    non_target_reads.append(pileupread.alignment)
                continue
            elif pileupread.indel > 0:
                indel_sequence = pileupread.alignment.query_sequence[
                    pileupread.query_position: pileupread.query_position + abs(pileupread.indel)
                ]
                if indel_sequence.upper() == target_base.upper():
                    target_reads.append(pileupread.alignment)
                else:
                    non_target_reads.append(pileupread.alignment)
                continue

            #Something else are SNPs
            #We have considered deletion and indel, so we filtered reads with None query_position
            if pileupread.query_position is not None:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                if base.upper() == target_base.upper():
                    target_reads.append(pileupread.alignment)
                else:
                    non_target_reads.append(pileupread.alignment)
    if invert == False:
        return target_reads
    else:
        return non_target_reads



def find_haplotype_specific_snps(bam_file_path, chrom, target_position, target_base):
    """
    The function here is to extend from the target position until none Y linked SNPs support and obtain all the Y linked reads.
    Directional extend up or down, separately by serving the target position as center.
    Get:
        1) ALL the Y linked reads
        2) The Y linked reads covered position range.
    """
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    
    target_reads = extract_reads_supporting_snp(bamfile, chrom, target_position-1, target_base)
    
    matched_reads = list(target_reads)
    unmatched_reads = list()
    
    min_position_up = extend_direction(bamfile, chrom, target_reads, target_position-1, "up", matched_reads, unmatched_reads)
    max_position_down = extend_direction(bamfile, chrom, target_reads, target_position-1, "down", matched_reads, unmatched_reads)
    
    match_out_bam="/Users/jiachangfu/Desktop/01_short_reads/03_Pret_aln/bam/test_matched.bam"
    unmatch_out_bam="/Users/jiachangfu/Desktop/01_short_reads/03_Pret_aln/bam/test_unmatched.bam"
    sorted_match_out_bam="/Users/jiachangfu/Desktop/01_short_reads/03_Pret_aln/bam/test_matched.sorted.bam"
    sorted_unmatch_out_bam="/Users/jiachangfu/Desktop/01_short_reads/03_Pret_aln/bam/test_unmatched.sorted.bam"

    with pysam.AlignmentFile(match_out_bam, "wb", header=bamfile.header) as out_bam:
        for read in matched_reads:
            out_bam.write(read)
    
    with pysam.AlignmentFile(unmatch_out_bam, "wb", header=bamfile.header) as out_bam:
        for read in unmatched_reads:
            out_bam.write(read)

    pysam.sort("-o", sorted_match_out_bam, match_out_bam)
    pysam.index(sorted_match_out_bam)

    pysam.sort("-o", sorted_unmatch_out_bam, unmatch_out_bam)
    pysam.index(sorted_unmatch_out_bam)
    

#    pysam.index(match_out_bam)
 #   pysam.index(unmatch_out_bam)
    
    bamfile.close()
    return matched_reads, min_position_up, max_position_down


def build_target_reference(bamfile, chrom, target_reads, quality_cutoff=0):
    """
    This function is aming at resolving the second round check of candidate informative SNPs.
    The first round filter is to make sure all the target reads have the same genotype, while the rest of reads in this region could also have shared genotype and non-shared genotype as well:
        This is because we need to extend the reads, so for non target reads shared the same variation could be served as potential extention reads.
        It works for the autosome don't  have any mutations between haploids. I mean, 0/0 or 1/1 for the autosomes, while the Y segements have the opposite genotype.
        But it failed to separate the situation that the female reference can be 0/1, and the Y genoetype is 0 or 1.
    So this second round check is to identify this false positive Y genotype by COMPARING THE FLANKING GENOTYPE of THESE READS, instead of only considerring one base.
    We firstly build a reference sequence based on target reference sequence. 
    For each reads could be served as candidate extension reads, we aligned them to the consensus reference to check the linked relationship.
    It breaks down into two steps. The first one is this function, build_target_reference, and the second is check_linked_mutations_with_reference.


    Details :
        1) We list all the target reads and obtain the upper and lower boundary bases.
        2) We only use high quality base to shape the common reference
        3) Bases with all the target reads support will be considered
    """
    
    min_pos = float("inf")
    max_pos = 0
    target_read_names = set(read.query_name for read in target_reads) 
    target_reads_seq = set(read.query_sequence for read in target_reads)
    for read in target_reads:
        min_pos = min(min_pos, read.reference_start)
        max_pos = max(max_pos, read.reference_end)
        
    #Target reference like this :  ('snp', 'T'), 5525608: ('snp', 'G'), 5525609: ('snp', 'T'), 5525610: ('snp', 'A'), 5525611: ('snp', 'A'), 5525612: ('snp', 'C'), 5525613: ('snp', 'T'), 5525614: ('snp', 'A'),   
    target_reference_mutations = {}
    
    for pileupcolumn in bamfile.pileup(chrom, start=min_pos, end=max_pos, min_base_quality=0, 
                                       min_mapping_quality=0, stepper="nofilter", truncate=True):
        position_mutations = {}
        
        for pileupread in pileupcolumn.pileups:
            
            #Only use target reads. The reads were selected both by reads name and reads sequence to avoid paired reads inconsistence
            if pileupread.alignment.query_name in target_read_names and pileupread.alignment.query_sequence in target_reads_seq:


                # Determine mutation type here
                if pileupread.is_del:
                    mutation_type = "deletion"
                    mutation_info = "-"
                elif pileupread.indel > 0:
                    mutation_type = "insertion"
                    mutation_info = pileupread.alignment.query_sequence[
                        pileupread.query_position : pileupread.query_position + pileupread.indel
                    ]
                elif pileupread.indel < 0:
                    mutation_type = "deletion"
                    mutation_info = "-"
                else:
                    mutation_type = "snp"
                    mutation_info = pileupread.alignment.query_sequence[pileupread.query_position]

                #for the reads have a deletion here, their pileupread.query_position could be None, we labelled them as - and pass them quality control        
                base_quality = 0
                if pileupread.query_position is None:
                    #print("None",pileupcolumn.reference_pos, pileupread.alignment.query_name,mutation_info, mutation_type)
                    base_quality == 60
                else:
                    base_quality = pileupread.alignment.query_qualities[pileupread.query_position]
                position_mutations[pileupread.alignment] = (mutation_type, mutation_info, base_quality)
        
        # Get consensus reference. We only use high quality base to construct targeted reads reference. And only keep the shared mutation among targeted reads to avoid potential wrong mapping or lower quality bases
        unique_mutations = set((mut[0], mut[1]) for mut in position_mutations.values() if mut[2] >= quality_cutoff )
        if len(unique_mutations) == 1:
            mutation_type, mutation_info = unique_mutations.pop()
            target_reference_mutations[pileupcolumn.reference_pos] = (mutation_type, mutation_info)
    return target_reference_mutations




def check_linked_mutations_with_reference(bamfile, chrom, target_reads, candidate_pos, target_reference, sum_check_threshold=2):
    """
    Here, we used reads supported candidate genotype as input , labelled as target read.
    Compare each reads to the reference constructed by the known high confident Y linked reads we identified before.
    We are tolerant to occasional sequencing error, so we used the sum check threshold 2 parameters.
    It can be tolerant with at most two reads with the same sequencing error. 
    If all the reads were Y supported, all these reads should share the same flanking mutations.
    Else they will be served as sharing the same genotype with one of the haplo group.
        


    Parameters
    ----------
    target_reads : TYPE
        They should be the reads shared the same candidate genotype. And here, by comparing the flanking region to the reference, we can determine its origin
    sum_check_threshold : See above description
        DESCRIPTION. The default is 2.

    """
    #Get target reads id and seq
    target_read_names = set(read.query_name for read in target_reads) 
    target_reads_seq = set(read.query_sequence for read in target_reads)

    #Only compare the shared region with reference
    ref_pos = list(target_reference.keys())
    for pileupcolumn in bamfile.pileup(chrom, start=min(target_reference.keys()), end =  max(target_reference.keys()) , min_base_quality=0, min_mapping_quality=0, stepper="nofilter", truncate=True):
        
        #Non-informative, just run over
        if pileupcolumn.reference_pos == candidate_pos:
            continue
        
        #Only focus on the shared region. While we constructed high quality reference, that means we have to filter some base with ambigous support. Some bases will be filter.
        if pileupcolumn.reference_pos not in ref_pos:
            continue
        
        sum_check = 0
        for pileupread in pileupcolumn.pileups:
            #Target reads filter
            if pileupread.alignment.query_name in target_read_names and pileupread.alignment.query_sequence in target_reads_seq:

                # Obtain the mutation type
                if pileupread.is_del:
                    mutation_type = "deletion"
                    mutation_info = "-"
                elif pileupread.indel > 0:
                    mutation_type = "insertion"
                    mutation_info = pileupread.alignment.query_sequence[
                        pileupread.query_position : pileupread.query_position + pileupread.indel
                    ]
                elif pileupread.indel < 0:
                    mutation_type = "deletion"
                    mutation_info = "-"
                else:
                    mutation_type = "snp"
                    mutation_info = pileupread.alignment.query_sequence[pileupread.query_position]
                
                #Ignore low quality bases
                quality = 0
                if pileupread.query_position is None:
                    quality = 60
                else:
                    quality = pileupread.alignment.query_qualities[pileupread.query_position]
                if (quality<20):
                    continue
                
                #Compare the reference and target reads
                ref_mutation = target_reference.get(pileupcolumn.reference_pos)
                #Debug
                #if pileupcolumn.reference_pos == 5525749 or  pileupcolumn.reference_pos == 5525921:
                    #print("ref",pileupcolumn.reference_pos,ref_mutation,mutation_type, mutation_info)
                
                if ref_mutation and (mutation_type, mutation_info) != ref_mutation:
                    sum_check+=1
        if sum_check >= sum_check_threshold:
            return False    
    return True

#def terminal_extend():
    #retrun()

def extend_direction(bamfile, chrom, target_reads, target_position, direction, matched_reads, unmatched_reads):
    """
    This is for extention the Y copy coordiate and get all the Y linked reads.
    It has down and up direction.
    It had two round check:
        1)Check the consistency between background bases and the targeted bases;
        2)Compare the flanking mutation of candidate mutation with the reference constructed by the high confident target reads;
    
    It truns back each reads and coordinate.
    """
    current_position = target_position

    while True:
        new_snps = []

        # target reads names to better intersect the pileup reads and our target reads
        target_read_names = set(read.query_name for read in target_reads)
        # since some paired reads have the same reads name ,so we added reads sequence as another identifier
        target_reads_seq = set(read.query_sequence for read in target_reads)

        #Up and Down direction check here
        start_pos = min(read.reference_start for read in target_reads) if direction == "up" else current_position  
        end_pos = current_position if direction == "up" else max(read.reference_end for read in target_reads) 

        print("\nstartPos: %s \n endpos : %s \n" %(start_pos, end_pos))
        for pileupcolumn in bamfile.pileup(chrom, start=start_pos, end=end_pos, truncate=True,
                                           stepper="nofilter",min_base_quality=0, min_mapping_quality=0):

            if (direction == "up" and pileupcolumn.pos < current_position) or \
               (direction == "down" and pileupcolumn.pos > current_position):

                bases = []
                non_target_bases = []
                for pileupread in pileupcolumn.pileups:                    
                    #Debug for some reads were empty match through below command                 if pileupread.alignment not in target_reads:
                    #if pileupcolumn.reference_pos <5526194 or pileupcolumn.reference_pos ==  5526291:
                       # print(pileupread.alignment)

                    #Targeted reads bases store to list
                    if pileupread.alignment.query_name in target_read_names and pileupread.alignment.query_sequence in target_reads_seq :
                        #Debug
#                        if pileupcolumn.reference_pos ==  5525921 or pileupcolumn.reference_pos ==  5526031 or pileupcolumn.reference_pos == 5525683:
#                            print("my aiming check for 5525683 : reads %s : is_del : %s. .indel = %s" % (pileupread.alignment.query_name , pileupread.is_del,pileupread.indel) )
                        
                        #Target reads 
                        if pileupread.is_del:
                            bases.append("-")
                        elif pileupread.indel > 0:
                            indel_seq = pileupread.alignment.query_sequence[
                                pileupread.query_position: pileupread.query_position + pileupread.indel
                                ]
                            if pileupcolumn.reference_pos ==  5525683:
                                print(indel_seq)
                        elif pileupread.indel < 0:
                            bases.append("-")

                        else:
                            base = pileupread.alignment.query_sequence[pileupread.query_position]
                            bases.append(base)
                            
                    #Non targeted reads bases store to list 
                    else:

                        if pileupread.is_del:
                            non_target_bases.append("-")
                        elif pileupread.indel > 0:
                            indel_seq = pileupread.alignment.query_sequence[
                                pileupread.query_position: pileupread.query_position + pileupread.indel
                                ]
                            non_target_bases.append(indel_seq)
                        elif pileupread.indel < 0:
                            non_target_bases.append("-")
                        else:
                            base = pileupread.alignment.query_sequence[pileupread.query_position]
                            non_target_bases.append(base)

                #list information to debug
                #print("\n",bases,name,pileupcolumn.reference_pos,temp)
                if pileupcolumn.reference_pos ==  5525683:
                    print("\nw\n",bases,non_target_bases)
                if len(bases) > 0 and bases.count(bases[0]) == len(bases): #in this 1bp region, the mutation of all associated targeted reads could be the same
                    target_genotype = bases[0]

                    #if pileupcolumn.reference_pos ==  5525683:
                      #  print("\nwaaa\n",bases,non_target_bases,name)

    #targeted bases ['T', 'T', 'T', 'T', 'T', 'T', 'T'] 
    #non targeted bases ['C', 'T', 'C', 'C', 'C', 'C', 'C', 'T', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'T', 'C', 'C'] 
    #Once there are some mutations differently in other reads, it could be served as 
                    ###The most important line : First round check to find candidate
                    ###The unmatched bases in non-targeted reads shoule be satisfied above 3 three times and there will be new reads can connect this area.
                    sum_check = sum(1 for non_target_base in non_target_bases if non_target_base not in bases  )
                    if sum_check > 2 and any(non_target_base in bases for non_target_base in non_target_bases) :
                            #print( "first_check : Position : %s target geno is %s ; non target geno is %s"  % (pileupcolumn.reference_pos, target_genotype, non_target_bases[0] ) )

                            target_reference = build_target_reference(bamfile, chrom, target_reads, quality_cutoff=20)
                            
                            test_target_reads = extract_reads_supporting_snp(bamfile, chrom, pileupcolumn.reference_pos, target_genotype)
                            if pileupcolumn.reference_pos ==  5525683:
                                print("\nwbbb\n",target_genotype, bases,non_target_bases,target_reference)
                            print(pileupcolumn.reference_pos )
                            #Second round check to examine the flanking mutation support
                            if check_linked_mutations_with_reference(bamfile, chrom, test_target_reads, pileupcolumn.reference_pos, target_reference):
                                new_snps.append((pileupcolumn.reference_pos, target_genotype))
                                print( "Position : %s target geno is %s ; non target geno is %s"  % (pileupcolumn.reference_pos, target_genotype, non_target_bases ) )
               # print("\ntargeted bases %s \n non targeted bases %s \n targeted_genotype : %s \n ref_position: %s " % (bases, non_target_bases, target_genotype,pileupcolumn.pos) )

        if not new_snps:
            break
        print(new_snps)
        new_target_position = min(new_snps, key=lambda x: x[0])[0] if direction == "up" else max(new_snps, key=lambda x: x[0])[0]
        target_genotype = min(new_snps, key=lambda x: x[0])[1] if direction == "up" else max(new_snps, key=lambda x: x[0])[1]
        print("\n new_target_postion %s \n new target genotype %s \n" % (new_target_position, target_genotype))

        target_reads = extract_reads_supporting_snp(bamfile, chrom, new_target_position, target_genotype)
        matched_reads.extend(target_reads)
        
        untarget_reads = extract_reads_supporting_snp(bamfile, chrom, new_target_position, target_genotype, invert=True)
        unmatched_reads.extend(untarget_reads)

        current_position = new_target_position

    return current_position


if __name__ == "__main__":
    bam_file_path = "/Users/jiachangfu/Desktop/01_short_reads/03_Pret_aln/bam/3.Pret_aln/NS.2118.001.IDT_i7_125---IDT_i5_125.114_R.q20.sorted.rmdup.bam"
    reference_fasta = "/Users/jiachangfu/Desktop/01_short_reads/0.reference_data/P.reticulata/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna"
    chrom = "NC_024336.1"
    #pos = 5524819 #position in vcf file, is 1 based
    target_base = "A"
    target_position = 5526293
   # target_position = 5526194
   # target_base = "G"

    
#    haplotypes = phase_reads(bam_path, reference_fasta, chrom, pos, target_base)
    final = find_haplotype_specific_snps(bam_file_path,chrom,target_position,target_base)      
    #for idx, hap in enumerate(haplotypes, 1):
     #   print(f"Haplotype {idx}: {hap}")





