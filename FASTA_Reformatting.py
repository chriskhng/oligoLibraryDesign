import re
import csv
import time
import shutil
from re import search
import numpy as np
import math

FASTA_file_name = "Phosphosite_PTM_seq.fasta"
entryInfo =[]
entryList = []
lstJoinedFASTASequences =[]
lstSequenceSeparated=[]
lstSplitEntry =[]
lst_masterFASTA = []
arb_num_FASTA = 0  # This is used as an on-off switch for this function(with a =+ at the end of the successful route)
while arb_num_FASTA == 0:
    try:
        start_time = time.time()  # Used later to report how long reformatting took.
        print("--- Reformatting FASTA entries ---")
        with open(FASTA_file_name) as FASTA_file:
            FASTA_data = FASTA_file.read()
        ## - File can close now
        lst_initialFASTA = FASTA_data.split('>GN:')
        # print(lst_initialFASTA)
        for entry in lst_initialFASTA[1:]:
            splitEntry =  entry.split('|')
            # print(splitEntry)
            lstSplitEntry.append(splitEntry)
            for item in lstSplitEntry:
                # print(item[3])
                sequenceSeparated= item[3].split('\n')
                # print(sequenceSeparated)
                geneID = sequenceSeparated.pop(0)
                sequenceSeparated.pop()
                joinedSequences =''.join(sequenceSeparated)
            entryInfo = splitEntry[0:3]
            # print(entryInfo)
            entryList = entryInfo
            entryList.append(geneID)
            entryList.append(joinedSequences)
            # print(entryList)
            lst_masterFASTA.append(entryList)
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        print("--- Reformatted FASTA entries in %s seconds ---" % (time.time() - start_time))
        arb_num_FASTA = 1
        print('First 3 entries of lst_masterFASTA:')
        print(lst_masterFASTA[0:3])
        print('number of entries:', len(lst_masterFASTA))
    except FileNotFoundError:
        FASTA_file_name = input("The FASTA file name is invalid, try again: ")
        ## - This list of sublists of protein entries will be used in Part 2.
        ## - This list will be searched for the user-input SoI to output the name of proteins containing the SoI

    FileName = 'CN reformatted PTM FASTA.csv'
    with open(FileName, 'w', newline='') as temp_FASTA_file:
        writer = csv.writer(temp_FASTA_file)
        writer.writerows(lst_masterFASTA)
    shutil.move(FileName, "generatedFiles/")