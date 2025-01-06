#!/usr/bin/env python

import os
import argparse as ap
import multiprocessing as mp

#Function for identifying the gap locations
def identify_gaps_in_proteome_alignment(contents):
    sequence_dict = {}
    key = ''

    #Variable for keeping track of the position of characters
    i = 0

    #List for storing the location of gaps that are found in any or all of the sequences
    gap_locations_list = []

    #Identify locations with gaps 
    for line in contents:
        #Store sequence name if it is encountered in sequence_dict as key
        if line.startswith('>'):
            key = line.strip()
            sequence_dict[key] = ''
            i = 0 #Reset position
        else:
            #Iterate through each character in a line
            for char in line.strip():
                sequence_dict[key] += char
                if char == '-' and i not in gap_locations_list:
                    gap_locations_list.append(i)
                i += 1 #Increment by 1 as we move from one character to the next

    #Convert dictionary to nested list to be used in the mapping latter 
    nested_seq_list = [[key, value] for key, value in sequence_dict.items()]
    return nested_seq_list, gap_locations_list

#Function for removing characters corresponding to gap locations
def remove_gaps(proteome_list):
    #Pull the proteome name 
    proteome_name = proteome_list[0]
    proteome_name = proteome_name[1:] #Simply remove the '>'

    print('Running the proteome:', proteome_name)

    #Pull the proteome sequence
    proteome_seq = proteome_list[1]

    gapless_proteome_seq = ''
    i = 0 
    
    #Removing the gaps 
    for char in proteome_seq:
        if i in gap_locations:
            pass 
        else:
            gapless_proteome_seq += char 
        i += 1

    return proteome_name, gapless_proteome_seq

if __name__ == '__main__':
    parser = ap.ArgumentParser(description='A script to remove gaps across the alignment')
    parser.add_argument('-i', '--input', help='Path to the input proteome alignment')
    parser.add_argument('-o', '--output', help='Path to write the output, gapless proteome alignment')
    args = parser.parse_args()

    #Open the proteome alignment file 
    alignment_path = args.input 
    with open(alignment_path, 'r') as input_alignment:
        alignment_contents = input_alignment.readlines()
    
    #Format the alignment into a nest list (for later multiprocessing) and identify the gap locations (in any and all sequences)
    nested_seq, gap_locations = identify_gaps_in_proteome_alignment(alignment_contents)

    #Run the removal of characters corresponding in gap locations (to any or all sequences) in parallel
    num_processes = 8 
    pool = mp.Pool(processes=num_processes)
    results = pool.map(remove_gaps, nested_seq)
    pool.close()
    pool.join()

    #Write the gapless sequence to single fasta 
    gapless_alignment_path = args.output
    with open(gapless_alignment_path, 'w') as outfile:
        for proteome_name, gapless_proteome_seq in results:
            outfile.write(f'>{proteome_name}\n{gapless_proteome_seq}\n')