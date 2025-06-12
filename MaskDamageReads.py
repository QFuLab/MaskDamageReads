# -*- coding: utf-8 -*-
#!/usr/bin/env python3

import pysam
import argparse

"""
MaskDamageReads.py

This script performs masking of potentially damaged bases in a BAM file 
either across the entire read or by a user-specified length from the 
ends, outputting a masked BAM file where damaged bases are changed to 
'N' and quality scores at these bases made 0. Masking is done in a 
strand aware manner, depending on whether the library is double or 
single strand. Please check the README for more information.

By: Dan JU & Fan BAI  

Citation: Fu Q, Cao P, Dai Q, Bennett EA, Feng X, Yang MA, et al. 
Denisovan mitochondrial DNA from dental calculus of the >146,000-year-old 
Harbin cranium. Cell. 2025
"""

def termini_par(value):
    try:
        return {'5p': int(value), '3p': int(value)}
    except ValueError:
        if value == 'all':
            return value
        elif ',' in value:
            (p5, p3) = value.split(',')
            return {'5p': int(p5), '3p': int(p3)}
        else:
            raise argparse.ArgumentTypeError(f"Invalid value: {value}. It must be an integer, 'all', or 'integer,integer'.")


def parse_args():
    parser = argparse.ArgumentParser(description="Masking the terminal damaged bases in BAM file. This script only works with bam file containing MD tag. The command `samtools calmd -b in.bam ref.fa > out.bam` can add MD tag for you.")
    
    parser.add_argument('-i', '--input', required=True, help="Input BAM file")
    parser.add_argument('-t', '--termini', type=termini_par, required=True, help="Number bp of the terminal bases to mask. If it is set to 'all', the script will mask the damaged bases in the whole read ('-t all' only works with '-m SS'). If the 5'' and 3'' ends need different base number, you need to seperate them by comma (no space between comma and the second number), and put 5'' number first. Should be attention that the 5'' and 3'' here refer to the reference direction (Value: an integer, 'all', or 'integer,integer')")
    parser.add_argument('-l', '--lib', choices=['SS', 'DS'], required=True, 
                        help="The way of library preparation: SS (single-strand library) or DS (double-strand library)")
    parser.add_argument('-o', '--output', required=True, help="Output BAM file")
    parser.add_argument('-m', '--mode', choices=['normal', 'low_coverage'], default='normal',
                        help="Mode selection: normal or low_coverage (Default: normal; see README for details)")
    parser.add_argument('--extract-damaged', action='store_true', 
                        help="Flag to extract damaged reads (Default: False)")
    parser.add_argument('--output-masked-position', action='store_true', 
                        help="Flag to output the coordinates of each N that this script added into a bed file (Default: False)")
    args = parser.parse_args()
    return args

def write_masked_position_to_bed(bam_out, masked_position_list):
    bed_fn = bam_out[:-4] + '.masked.position.bed'
    with open(bed_fn, 'w') as f:
        for ref_id, ref_pos in masked_position_list:
            f.write(f'{ref_id}\t{ref_pos}\t{int(ref_pos) + 1}\n')

def mask_base_in_reads(read, read_seq, read_qual, lib, termini, qpos, read_base, ref_base):
    read_len = len(read_seq)
    base_is_masked = False

    if lib == 'SS':
        if read.is_forward and (qpos < termini['5p'] or read_len - qpos <= termini['3p']):
            # Mask 5'' and 3'' terminal C to T in SS masking
            if ref_base == 'C' and read_base == 'T':
                read_seq[qpos] = 'N'
                if read_qual is not None:
                    read_qual[qpos] = 0  
                base_is_masked = True  
        elif read.is_reverse and (qpos < termini['3p'] or read_len - qpos <= termini['5p']):
            if ref_base == 'G' and read_base == 'A':
                read_seq[qpos] = 'N' 
                if read_qual is not None:
                    read_qual[qpos] = 0  
                base_is_masked = True  
    elif lib == 'DS':
        # The qpos offsets from start of reads. For forward strand, it starts from 5'',
        # but for reverse strand, it starts from 3''. However, the reverse strand is stored in
        # complementary based.
        if read.is_forward:
            if qpos < termini['5p'] and ref_base == 'C' and read_base == 'T':
                # This is 5'' C to T of forward strand, and also 3'' G to A of reverse strand
                read_seq[qpos] = 'N'  
                if read_qual is not None:
                    read_qual[qpos] = 0   
                base_is_masked = True    
            if read_len - qpos <= termini['3p'] and ref_base == 'G' and read_base == 'A':
                # Same as above one
                read_seq[qpos] = 'N' 
                if read_qual is not None:
                    read_qual[qpos] = 0  
                base_is_masked = True     
        elif read.is_reverse:
            if qpos < termini['3p'] and ref_base == 'C' and read_base == 'T':
                # This is 5'' C to T of forward strand, and also 3'' G to A of reverse strand
                read_seq[qpos] = 'N'  
                if read_qual is not None:
                    read_qual[qpos] = 0   
                base_is_masked = True    
            if read_len - qpos <= termini['5p'] and ref_base == 'G' and read_base == 'A':
                # Same as above one
                read_seq[qpos] = 'N' 
                if read_qual is not None:
                    read_qual[qpos] = 0  
                base_is_masked = True     
    return (read_seq, read_qual, base_is_masked)

def mask_damaged(input_bam, output_bam, termini, lib, output_masked_position):
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, 'wb', template=bam_in)

    masked_position_list = []
    # Iterate over each read in BAM
    for read in bam_in.fetch():
        # If unmapped write as is
        if read.is_unmapped:
            bam_out.write(read)
            continue
            
        read_seq = list(read.query_sequence)
        read_qual = read.query_qualities

        # Reference free version
        # Get the aligned pairs of reference and query positions
        # The reference bases are got from MD tag
        for qpos, refpos, refbase in read.get_aligned_pairs(matches_only=True, with_seq=True):
            if refpos is not None and qpos is not None:
                read_base = read_seq[qpos].upper()
                ref_base = refbase.upper()
                read_seq, read_qual, base_is_masked = mask_base_in_reads(read, read_seq, read_qual, lib, termini, qpos, read_base, ref_base)
                if base_is_masked:
                    masked_position_list.append((bam_in.get_reference_name(read.reference_id), refpos))

        # Update and write to output BAM masked sequence
        read.query_sequence = ''.join(read_seq)
        # Add back original quality scores
        read.query_qualities = read_qual
        
        bam_out.write(read)

    if output_masked_position:
        write_masked_position_to_bed(output_bam, masked_position_list)

    # Close input and output BAM files
    bam_in.close()
    bam_out.close()

    print('Finish')
    return 0

def SS_lib(input_bam, output_bam, output_masked_position):
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, 'wb', template=bam_in)

    masked_position_list = []
    # Iterate over each read in BAM
    for read in bam_in.fetch():
        # If unmapped write as is
        if read.is_unmapped:
            bam_out.write(read)
            continue
            
        read_seq = list(read.query_sequence)
        read_qual = read.query_qualities

        # Reference free version
        # Get the aligned pairs of reference and query positions
        # The reference bases are got from MD tag
        for qpos, refpos, refbase in read.get_aligned_pairs(matches_only=True, with_seq=True):
            if refpos is not None and qpos is not None:
                read_base = read_seq[qpos].upper()
                ref_base = refbase.upper()
                # Mask 5'' and 3'' terminal C to T in SS lib
                if read.is_forward and ref_base == 'C' and read_base == 'T':
                    read_seq[qpos] = 'N'  
                    if read_qual is not None:
                        read_qual[qpos] = 0 
                    masked_position_list.append((bam_in.get_reference_name(read.reference_id), refpos))
                elif read.is_reverse and ref_base == 'G' and read_base == 'A':
                    read_seq[qpos] = 'N' 
                    if read_qual is not None:
                        read_qual[qpos] = 0 
                    masked_position_list.append((bam_in.get_reference_name(read.reference_id), refpos))

        # Update and write to output BAM masked sequence
        read.query_sequence = ''.join(read_seq)
        # Add back original quality scores
        read.query_qualities = read_qual
        
        bam_out.write(read)

    if output_masked_position:
        write_masked_position_to_bed(output_bam, masked_position_list)

    # Close input and output BAM files
    bam_in.close()
    bam_out.close()

    print('Finish')
    return 0

def mask_base_in_reads_lc(read, read_seq, read_qual, lib, termini, qpos, read_base, ref_base):
    read_len = len(read_seq)
    base_is_masked = False

    if lib == 'SS':
        if read.is_forward and (qpos < termini['5p'] or read_len - qpos <= termini['3p']):
            # Mask 5'' and 3'' terminal C to T in SS masking
            if (read_base == 'T' or read_base == 'C'):
                read_seq[qpos] = 'N'
                if read_qual is not None:
                    read_qual[qpos] = 0  
                base_is_masked = True  
        elif read.is_reverse and (qpos < termini['3p'] or read_len - qpos <= termini['5p']):
            if (read_base == 'A' or read_base == 'G'):
                read_seq[qpos] = 'N' 
                if read_qual is not None:
                    read_qual[qpos] = 0  
                base_is_masked = True  
    elif lib == 'DS':
        # The qpos offsets from start of reads. For forward strand, it starts from 5'',
        # but for reverse strand, it starts from 3''. However, the reverse strand is stored in
        # complementary based.
        if read.is_forward:
            if qpos < termini['5p'] and (read_base == 'T' or read_base == 'C'):
                # This is 5'' C to T of forward strand, and also 3'' G to A of reverse strand
                read_seq[qpos] = 'N'  
                if read_qual is not None:
                    read_qual[qpos] = 0   
                base_is_masked = True    
            if read_len - qpos <= termini['3p'] and (read_base == 'A' or read_base == 'G'):
                # Same as above one
                read_seq[qpos] = 'N' 
                if read_qual is not None:
                    read_qual[qpos] = 0  
                base_is_masked = True     
        elif read.is_reverse:
            if qpos < termini['3p'] and (read_base == 'T' or read_base == 'C'):
                # This is 5'' C to T of forward strand, and also 3'' G to A of reverse strand
                read_seq[qpos] = 'N'  
                if read_qual is not None:
                    read_qual[qpos] = 0   
                base_is_masked = True    
            if read_len - qpos <= termini['5p'] and (read_base == 'A' or read_base == 'G'):
                # Same as above one
                read_seq[qpos] = 'N' 
                if read_qual is not None:
                    read_qual[qpos] = 0  
                base_is_masked = True     
    return (read_seq, read_qual, base_is_masked)

def mask_damaged_lc(input_bam, output_bam, termini, lib, output_masked_position):
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, 'wb', template=bam_in)

    masked_position_list = []
    # Iterate over each read in BAM
    for read in bam_in.fetch():
        # If unmapped write as is
        if read.is_unmapped:
            bam_out.write(read)
            continue
            
        read_seq = list(read.query_sequence)
        read_qual = read.query_qualities

        # Reference free version
        # Get the aligned pairs of reference and query positions
        # The reference bases are got from MD tag
        for qpos, refpos, refbase in read.get_aligned_pairs(matches_only=True, with_seq=True):
            if refpos is not None and qpos is not None:
                read_base = read_seq[qpos].upper()
                ref_base = refbase.upper()
                read_seq, read_qual, base_is_masked = mask_base_in_reads_lc(read, read_seq, read_qual, lib, termini, qpos, read_base, ref_base)
                if base_is_masked:
                    masked_position_list.append((bam_in.get_reference_name(read.reference_id), refpos))

        # Update and write to output BAM masked sequence
        read.query_sequence = ''.join(read_seq)
        # Add back original quality scores
        read.query_qualities = read_qual
        
        bam_out.write(read)

    if output_masked_position:
        write_masked_position_to_bed(output_bam, masked_position_list)

    # Close input and output BAM files
    bam_in.close()
    bam_out.close()

    print('Finish')
    return 0

def SS_lib_lc(input_bam, output_bam, output_masked_position):
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, 'wb', template=bam_in)

    masked_position_list = []
    # Iterate over each read in BAM
    for read in bam_in.fetch():
        # If unmapped write as is
        if read.is_unmapped:
            bam_out.write(read)
            continue
            
        read_seq = list(read.query_sequence)
        read_qual = read.query_qualities

        # Reference free version
        # Get the aligned pairs of reference and query positions
        # The reference bases are got from MD tag
        for qpos, refpos, refbase in read.get_aligned_pairs(matches_only=True, with_seq=True):
            if refpos is not None and qpos is not None:
                read_base = read_seq[qpos].upper()
                ref_base = refbase.upper()
                # Mask 5'' and 3'' terminal C to T in SS lib
                if read.is_forward and (read_base == 'T' or read_base == 'C'):
                    read_seq[qpos] = 'N'  
                    if read_qual is not None:
                        read_qual[qpos] = 0 
                    masked_position_list.append((bam_in.get_reference_name(read.reference_id), refpos))
                elif read.is_reverse and (read_base == 'A' or read_base == 'G'):
                    read_seq[qpos] = 'N' 
                    if read_qual is not None:
                        read_qual[qpos] = 0 
                    masked_position_list.append((bam_in.get_reference_name(read.reference_id), refpos))

        # Update and write to output BAM masked sequence
        read.query_sequence = ''.join(read_seq)
        # Add back original quality scores
        read.query_qualities = read_qual
        
        bam_out.write(read)

    if output_masked_position:
        write_masked_position_to_bed(output_bam, masked_position_list)

    # Close input and output BAM files
    bam_in.close()
    bam_out.close()

    print('Finish')
    return 0

def terminal_base_is_damaged(read, read_seq, lib, termini, qpos, read_base, ref_base):
    read_len = len(read_seq)

    if lib == 'SS':
        if read.is_forward:
            if qpos < termini['5p'] or read_len - qpos <= termini['3p']:
                # 5'' and 3'' terminal C to T in SS lib is damaged base
                if ref_base == 'C' and read_base == 'T':
                    return True 
        elif read.is_reverse:
            if qpos < termini['3p'] or read_len - qpos <= termini['5p']:
                if ref_base == 'G' and read_base == 'A':
                    return True
    elif lib == 'DS':
        # The qpos offsets from start of reads. For forward strand, it starts from 5'',
        # but for reverse strand, it starts from 3''. However, the reverse strand is stored in
        # complementary based.
        if read.is_forward:
            if qpos < termini['5p'] and ref_base == 'C' and read_base == 'T':
                # This is 5'' C to T of forward strand, and also 3'' G to A of reverse strand
                return True
            if read_len - qpos <= termini['3p'] and ref_base == 'G' and read_base == 'A':
                return True             
        elif read.is_reverse:
            if qpos < termini['3p'] and ref_base == 'C' and read_base == 'T':
                return True             
            if read_len - qpos <= termini['5p'] and ref_base == 'G' and read_base == 'A':
                return True     
    
    return False

def read_is_damaged(read, termini, lib):
    if read.is_unmapped:
        return False
      
    read_seq = list(read.query_sequence)

    # Reference free version
    # Get the aligned pairs of reference and query positions
    # The reference bases are got from MD tag
    for qpos, refpos, refbase in read.get_aligned_pairs(matches_only=True, with_seq=True):
        if refpos is not None and qpos is not None:
            read_base = read_seq[qpos].upper()
            ref_base = refbase.upper()
            if terminal_base_is_damaged(read, read_seq, lib, termini, qpos, read_base, ref_base):
                return True
    
    return False

def extract_damaged(input_bam, output_bam, termini, lib):
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, 'wb', template=bam_in)

    # Iterate over each read in BAM
    for read in bam_in.fetch():
        # If damaged write as is
        if read_is_damaged(read, termini, lib):
            bam_out.write(read)

    # Close input and output BAM files
    bam_in.close()
    bam_out.close()

    print('Finish')
    return 0

def main():
    args = parse_args()
    
    print('==============PARAMETER===================')
    print(f'Input BAM: {args.input}')
    print(f'Output BAM: {args.output}')
    print(f'Termini (bp to mask): {args.termini}')
    print(f'Library type: {args.lib}')
    print(f'Masking mode: {args.mode}')
    print(f'Extract damage: {args.extract_damaged}')
    print(f'Output masked position: {args.output_masked_position}')
    print('==============START===================')

    if args.extract_damaged is False and args.mode == 'normal':
        print(f'Mask {args.input} output to {args.output} ... ', sep=' ')
        if args.termini != 'all':
            mask_damaged(args.input, args.output, args.termini, args.lib, args.output_masked_position)
        elif args.lib == 'SS' and args.termini == 'all':
            SS_lib(args.input, args.output, args.output_masked_position)
        else:
            print("ERROR: '-t all' only works with '-m SS'")
            return 1
    elif args.extract_damaged is False and args.mode == 'low_coverage':
        print(f'Mask {args.input} output to {args.output} ... ', sep=' ')
        if args.termini != 'all':
            mask_damaged_lc(args.input, args.output, args.termini, args.lib, args.output_masked_position)
        elif args.lib == 'SS' and args.termini == 'all':
            SS_lib_lc(args.input, args.output, args.output_masked_position)
        else:
            print("ERROR: '-t all' only works with '-m SS'")
            return 1
    else:
        print(f'Extract damaged reads from {args.input} output to {args.output} ... ', sep=' ')
        extract_damaged(args.input, args.output, args.termini, args.lib)

    return 0

if __name__ == '__main__':
    main()