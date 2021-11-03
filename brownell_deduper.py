#!/usr/bin/env python3

import argparse
import re


#argparse 
#randomers are used automatically if an umi file is not provided
#paired end files can be taken if --paired is declared
#duplicate choosing by total qscore can be enabled with -q
#outputfile can be named, but defaults 

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="path to file")
parser.add_argument("-p", "--paired", default = False, action = 'store_true', help = "Declare this option to indicate paired-end data.")
parser.add_argument("-u", "--umi", default = False, help = "Designate location of umi file. Defaults to randomers.")
parser.add_argument("-o", "--output_filename", default = False, help = "Give an optional filename for the deduplicated output. Defaults to adding '_deduped' to the input filename.")
parser.add_argument("-q", "--select_highest_total_quality", default = False, action = 'store_true', help = "Use this option to select highest total quality duplicate. False by default.")


args = parser.parse_args()

filepath = args.file
paired = args.paired
umifile = args.umi
OUTPUT_NAME = args.output_filename
select_highest_sum = args.select_highest_total_quality





#outputfile, inputfile and umifile are opened. Umifile is parsed.   
if not OUTPUT_NAME:
    outputfile = open(filepath.strip(".sam") + "_deduped.sam", "w")
else:
    outputfile = open(OUTPUT_NAME, "w") 

if umifile:
    umiset = set()
    with open(umifile, "r") as known_umis:
        for line in known_umis:
            umiset.add(line.strip())

sam = open(filepath, "r")








#below are the three functions used in the main loop. Most are self-explanatory

def find_umi(input):
    return input.split(":")[-1]

#adjust position is the function used to parse the cigar string
#if the input strand is +, then it simply subtracts soft-clipping on the left fron the position
#if the input strand is -, then it adds up all numbers in the string, and subtracts the number of indexes and soft clipping on the left. It also subtracts one. 
#This process is intended to get the actual adjusted rightmost mapping position for - strand reads. 

def adjust_position(strand_is_minus, position, cigar):
    cig_split = re.split("(\D)", cigar)
    increment = 0
    if strand_is_minus:
        just_digits = tuple(map(float, filter(str.isdigit, cig_split)))
        increment += sum(just_digits) - 1
        if cig_split[1] == "S" or cig_split[1] == "H":
            increment -= just_digits[0]
        if "I" in cig_split:
            for i, element in enumerate(cig_split):
                if element == "I":
                    increment -= float(cig_split[i-1])
    else:
        if cig_split[1] == "S" or cig_split[1] == "H":
            increment -= int(cig_split[0])
    return position + increment 
        





def minus_strand(flag):
    return (int(flag) & 16) == 16
        
        
def is_pair_1(flag):
    return (int(flag) & 64) == 64 


def total_quality(string):
    return sum(ord(x) - 33 for x in string.strip())
    



#The below variables are used in the main loop to store values. 
#Last_seen_chromosome stores whatever the newest chromosome is.
#uniques stores whatever unique lines are within the current chromosome. When a new chromosome is enocountered, the lines are dumped to file. 

total = 0
headers = 0
grand_total = 0
last_seen_chromosome = ""
dupcounter = 0
total_uniques = 0
unknown_umis = 0

qscoredict = {}
uniques = {}
pairbank = {}

for line in sam:
    
    grand_total += 1
    if line[0] != '@':
        total += 1
        
        #the lines are each split by tab, then their important features are assigned to variables
        #these important features are umi, chromosome, bitwise real start position, etc.

        linesep = line.split("\t")
        qname = linesep[0]
        umi = find_umi(qname)

        if umifile and umi not in umiset:
            #skipping over bad umis
            unknown_umis += 1
            continue 
        


        

        
        chromosome = linesep[2]
        bitwise = linesep[1]
        
        
        is_minus_strand = minus_strand(linesep[1])
        cigar = linesep[5]

        if paired:
            pair1 = is_pair_1(bitwise)
        
        #the position is adjusted for soft clipping and other position-altering elements of the alignment
        position = adjust_position(is_minus_strand, int(linesep[3]), cigar)

        #Below is where files are written to output if a new chromosome is encountered. 
        #Last seen chromosome is updated to the new one, if so. 
        #Comparing within chromosomes limits the amount of memory in use at a given time without affecting comparisons. 
        if chromosome != last_seen_chromosome:
            outputfile.write("".join(uniques.values()))
            print("finished chromosome " + last_seen_chromosome + ". Uniques: " + str(len(uniques)))
            total_uniques += len(uniques)
            last_seen_chromosome = chromosome

            #certain dictionaries associated with other options are reset as well
            uniques = {}
            pairbank = {}
            qscoredict = {}
        
        #seqinfo is the key we use to determine uniqueness. It comtains all relevant information for determining PCR duplicates, except chromosome. 
        #this is because chromosome is already the same, due to the input file being sorted and the dumping behavior. 
        seqinfo = (umi, is_minus_strand, position)

        #If not paired, then the logic below skips over duplicates already logged in the uniques dictionary
        #...or it logs a newly seen read into said dictionary 

        #Note that the lines of #### denote the nestedness of the quality selector option. \
        if not paired:    
            if seqinfo in uniques:
                dupcounter += 1

                #############################################################
                if select_highest_sum:
                    current_line_sum_quality = total_quality(linesep[10])
                    if current_line_sum_quality > qscoredict[seqinfo]:
                        uniques[seqinfo] = line
                        qscoredict[seqinfo] = current_line_sum_quality
                
            


            else:  
                uniques[seqinfo] = line
                if select_highest_sum:
                    qscoredict[seqinfo] = total_quality(linesep[10])

        else:
            #If the input is paired, an alternative logic is used. 
            #Lines of a mate-pair (defined by matching qname) are logged into the pairbank (a dictionary where single pairs stay until the mate is encountered).
            #Once the other line of the mate-pair is encountered, the two lines are checked against the uniques dictionary using both sets of information
            #As with before, mate-pairs that already have their seqinfo (umi, is_minus_strand, position) already logged are discarded
            #Otherwise, they are logged. 
            

            if qname in pairbank:
                
                if pair1:
                    
                    if (seqinfo, pairbank[qname][0]) in uniques:
                        dupcounter += 2
                        #####################################################################################
                        if select_highest_sum:
                            pair_total_quality = pairbank[qname][2] + total_quality(linesep[10])
                            if pair_total_quality > qscoredict[(seqinfo, pairbank[qname][0])]:
                                uniques[(seqinfo, pairbank[qname][0])] = ("".join((line, pairbank[qname][1])))
                                qscoredict[(seqinfo, pairbank[qname][0])] = pair_total_quality
                    else:
                        uniques[(seqinfo, pairbank[qname][0])] = ("".join((line, pairbank[qname][1])))
                        if select_highest_sum:
                            qscoredict[(seqinfo, pairbank[qname][0])] = pairbank[qname][2] + total_quality(linesep[10])
                else:
                    
                    if (pairbank[qname][0], seqinfo) in uniques:
                        dupcounter += 2
                        #####################################################################################
                        if select_highest_sum:
                            pair_total_quality = pairbank[qname][2] + total_quality(linesep[10])
                            if pair_total_quality > qscoredict[(pairbank[qname][0], seqinfo)]:
                                uniques[(pairbank[qname][0], seqinfo)] = ("".join((pairbank[qname][1], line)))
                                qscoredict[(pairbank[qname][0], seqinfo)] = pair_total_quality
                        
                    else:
                        uniques[(pairbank[qname][0], seqinfo)] = ("".join((pairbank[qname][1], line)))
                        if select_highest_sum:
                            
                            qscoredict[(pairbank[qname][0], seqinfo)] = pairbank[qname][2] + total_quality(linesep[10])

            else:
                
                pairbank[qname] = (seqinfo, line) if not select_highest_sum else (seqinfo, line, total_quality(linesep[10]))
                
            
    else:
        #headers are written to the file. 
        outputfile.write(line)
        headers += 1
            
#remaining unique lines left in the uniques dictionary are written to file
outputfile.write("".join(uniques.values()))
total_uniques += len(uniques)
print("finished chromosome " + last_seen_chromosome + ". Uniques: " + str(len(uniques)))




#Some stats are printed

print("unknown umis")
print(unknown_umis)

print("total lines:")
print(grand_total)

print("total headers:")
print(headers)

print("duplicates:")          
print(dupcounter)

print("total uniques:")
print(total_uniques)


        

