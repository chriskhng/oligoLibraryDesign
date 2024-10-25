import re
import csv
import time
import shutil
from re import search

import numpy as np
import math

# FASTA_file_name = "Phosphosite_PTM_seq.fasta"

#
#
# def FASTA_reformat(FASTA_file_name):
#     lst_masterFASTA = []
#     arb_num_FASTA = 0 # This is used as an on-off switch for this function (there is a =+ at the end of the successful route)
#     while arb_num_FASTA == 0:
#         try:
#             start_time = time.time() #Used later to report how long reformatting took.
#             print("--- Reformatting FASTA entries ---")
#
#             with open(FASTA_file_name) as FASTA_file:
#                 FASTA_data = FASTA_file.read()
#             ## - File can close now
#             FASTA_newlineM = FASTA_data.replace('\nM', '|M') #this replaces the '/n' delimits before the sequence to match the rest of the '|' delimits.
#             FASTA_newlineReplaced = FASTA_newlineM.replace('\n','') #this removes the '/n' within the sequences
#             FASTA_commaDelimited = FASTA_newlineReplaced.replace('|',',') #this replaces |  with ,
#             ## - lst_masterFASTAwithAuthors is a list of each protein entry, containing protein name, sequence, etc, in a single string.
#             lst_masterFASTAwithAuthors = FASTA_commaDelimited.split('>GN:')   ## I split the starting FASTA with \n>GN as that's what's separtating the each entries in the FASTA:
#             ## - lst_masterFASTAwithAuthors is a list of each protein entry: a list containing protein name, sequence, etc, in a single string.
#             lst_rawFASTA = lst_masterFASTAwithAuthors[1:] ##This removes the author watermark at the top which is index 0 here.
#             ## -Below, relevant info will be pulled out from each protein entry
#             ## - We are only interested in the gene name and protein sequence
#             for protein_entry in lst_rawFASTA:
#             ## - The gene name is pulled from each FASTA entry and saved into FASTgeneNames
#             ## - This is done by indexing from 0 to the first comma (which delimits name from the next entry):
#                 FASTAgeneNames = protein_entry[0: protein_entry.find(',')]
#                 proteinSeqFinal = protein_entry[protein_entry.find(',M')+1: ] #This is pulling the sequence into proteinSeqFinal by indexing from Methionine to the end
#
#             ## - Now, we put protein name and protein sequence together into a list for each entry
#             ## - First, the protein names are separated into genomic name, common name, and ID name
#                 lst_FASTAentry = FASTAgeneNames.split(' ')
#             ## - Next, the protein sequence is appended to end of the names.
#                 lst_FASTAentry.append (proteinSeqFinal)
#             ## - Finally, the list for the entry is appended to the growing list of protein entries from the FASTA.
#                 lst_masterFASTA.append (lst_FASTAentry)
#
#
#             print("--- Reformatted FASTA entries in %s seconds ---" % (time.time() - start_time))
#             arb_num_FASTA = 1
#             print('First 3 entries of lst_masterFASTA:')
#             print(lst_masterFASTA[0:3])
#
#         except FileNotFoundError:
#             FASTA_file_name = input ("The FASTA file name is invalid, try again: ")
#             ## - This list of sublists of protein entries will be used in Part 2.
#             ## - This list will be searched for the user-input SoI to output the name of proteins containing the SoI
#     return lst_masterFASTA

# FASTA = FASTA_reformat("Phosphosite_PTM_seq.fasta")
# print(FASTA[2:4])



reformatted_FASTA_file_name = "generatedFiles/CN reformatted PTM FASTA.csv"
def openingReformattedFASTA(reformatted_FASTA_file_name):
    with open(reformatted_FASTA_file_name) as FASTA_file:
        FASTA_data = FASTA_file.read()
    # print(FASTA_data[:10000])
    delimited_by_newLine = FASTA_data.split('\n') #Spliting the csv into entry by /n. each entry is one string
    reformatted_FASTA_List = []
    for entry in delimited_by_newLine: #This splits each entry (a string) by the ',' and puts it in a growing lst
        entry.split(',')
        reformatted_FASTA_List.append(entry.split(','))
    # print(delimited_by_newLine[:5])
    print(reformatted_FASTA_List[:10])
    return reformatted_FASTA_List

# user_input_SoI = "s..s..s"
lst_SoI_results = []
# pattern = 's..s..s'
# inputPlusMinus = 15
# spanPlusMinus = inputPlusMinus - math.ceil(len(pattern)/2)+1

def plusMinusSoI_search_FASTA(FASTA,pattern,inputPlusMinus):
    spanPlusMinus = inputPlusMinus - math.ceil(len(pattern) / 2) + 1 # this calculates how many nucleotides flank the pattern of interest for a given plus-minus
    counterSoI = 0
    counterProteins = 0
    lst_SoI_results = []

    for protein_entry in FASTA[:-1]: # it is :-1 because the reformated FASTA has an empty lst at the end. so I did this sphagetti code.
        # print(protein_entry[0:5])
        searchSpan = re.search(pattern, protein_entry[4]) #this searchs to see if the FASTA sequence has a match to the input pattern
        # print(searchSpan)
        if searchSpan: #if search returns something
            matches = re.findall(pattern, protein_entry[4]) #This checks to see how many matches there are in the sequence. matches is a list of the hit(s)
            if len(matches) > 1:  #if there is more than one hit
                # print('more than one')
                # print(matches)
                for multiPattern in matches: #for each hit in all the hits in the sequence
                    multiSearchSpan = re.search(multiPattern, protein_entry[4]) #This finds the indices of the hit
                    # print(multiSearchSpan)
                    # print(multiPattern)
                    #multiLowRange and HighRange are the indices of the peptide induceing the plus minus
                    #the np.clip accounts for if the
                    multiLowRange = np.clip(multiSearchSpan.start() - (spanPlusMinus), 0, None)
                    multiHighRange = multiSearchSpan.end() + (spanPlusMinus)
                    lst_SoI_results.append([protein_entry[0], protein_entry[4][
                                                              multiLowRange:multiHighRange],multiPattern,protein_entry[1],protein_entry[2], protein_entry[3]])
                    counterSoI += 1
            else:
                lowRange = np.clip(searchSpan.start() - (spanPlusMinus), 0, None)
                highRange = searchSpan.end() + (spanPlusMinus)
                sequenceHit = re.search(pattern, protein_entry[4])
                # print(sequenceHit.group())
                # print(protein_entry[1][lowRange:highRange])  ##THIS IS THE ONE
                lst_SoI_results.append([protein_entry[0],protein_entry[4][
                                                          lowRange:highRange],sequenceHit.group(),protein_entry[1],protein_entry[2], protein_entry[3]])
                counterSoI += 1
            counterProteins +=1
        else:
            continue
    print('***** Sequence of Interest search completed for', pattern,'for plus-minus', inputPlusMinus,'*****')
    print('Number of Proteins:', counterProteins)
    print('Number of SoI detected:',counterSoI)

    return lst_SoI_results

def sToJSwap(plusMinusOutput,pattern):
    lst_s2J_results = []
    for entry in plusMinusOutput:
        # print("entry:", entry)
        # print(re.search(entry[2], entry[1]))  ## This searchs for the sequence hit in the plus mins sequence
        firstSerineIndex = (re.search(entry[2], entry[1]).start())  ### first serine
        # print('firstSerineIndex:', firstSerineIndex)
        thirdSerineIndex = (re.search(entry[2], entry[1]).end()) - 1  #### Third serine
        # print('thirdSerineIndex:', thirdSerineIndex)
        ###central serine
        # print('re.research output with [1:]:', re.search('s', pattern[1:])) # re.findall output: <re.Match object; span=(2, 3), match='s'>
        centralSerineIndexTemp = re.search('s', pattern[1:])
        # print('centralSerineIndexTemp:', centralSerineIndexTemp)
        # print('centralSerineIndexTemp.start():', centralSerineIndexTemp.start())
        centralSerineIndexAdd = centralSerineIndexTemp.start() + 1
        # print('centralSerineAdd:', centralSerineIndexAdd)
        centralSerineIndex = centralSerineIndexAdd + firstSerineIndex
        # print('centralSerineIndex:', centralSerineIndex)

        seqIntoList = list(entry[1])
        # print(seqIntoList)
        seqIntoList[firstSerineIndex] = ("J")
        seqIntoList[centralSerineIndex] = ("J")
        seqIntoList[thirdSerineIndex] = ("J")
        # print(seqIntoList)
        # print(tempList[17])
        rejoinedList = ''.join(seqIntoList)
        # print(rejoinedList)
        # print('')
        lst_s2J_results.append([entry[0], rejoinedList, entry[2], entry[3], entry[4],entry[5]])
    return lst_s2J_results

