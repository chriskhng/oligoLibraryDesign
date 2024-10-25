import re
import csv
import time
import shutil
import numpy as np
import math

from oligoLibraryDesignFunctions import plusMinusSoI_search_FASTA, sToJSwap, openingReformattedFASTA

time_str = time.strftime("%Y%m%d-%H%M%S") ## - Used to prefix output files.

### Inputs ###
inputPattern = 'sss'
inputPlusMinus = 15
reformatted_FASTA_file_name = "generatedFiles/CN reformatted PTM FASTA.csv"

startingNucleotides = 'TTATAATCATCCTCCCCGGCGGTACCAAG'
endingNucleotides = 'AAGCTTTAACGAGCACACATCCTATTTGG'

restrictionSite01 = 'ACCAAG' ###DpnI
restrictionSite02 = 'AAGCTT' ###HindIII
restrictionKSF = 'AAAAGCTTT'
restrictionESF = 'GAAAGCTTT'
fixedKSF = 'AAAAGTTTT'
fixedESF = 'GAAAGTTTT'


# FASTA = FASTA_reformat("Phosphosite_PTM_seq.fasta")
# print('')
FASTA = openingReformattedFASTA(reformatted_FASTA_file_name)



spanPlusMinus = inputPlusMinus - math.ceil(len(inputPattern)/2)+1
plusMinusOutput = plusMinusSoI_search_FASTA(FASTA,inputPattern,inputPlusMinus)
print("With", inputPattern, "and PlusMinus", inputPlusMinus)
print(plusMinusOutput[0:20])
print('Number of entries found:' ,len(plusMinusOutput))


### s to J swaps ###
print('')
lst_s2J_results = sToJSwap(plusMinusOutput,inputPattern)
print('--- s to J swaps done --- ')
print('here are the first few:')
print(lst_s2J_results[0:20])
print('number of items:', len(lst_s2J_results))
print('')


codonDictionary = {
    'A':'GCG','R':'CGC','N':'AAC','D': 'GAT','C': 'TGC',
    'E': 'GAA','Q': 'CAG','G': 'GGC','H': 'CAT','I':'ATT',
    'L':'CTG','K': 'AAA','M': 'ATG','F': 'TTT','P': 'CCG',
    'S': 'AGC','T': 'ACC','W': 'TGG','Y': 'TAT','V': 'GTG',
    's': 'AGC','J': 'TAG'
}


# aaSequenceSample = 'ttPsEsPRAQAtJRLJtAJcPtPKVQsRCss'
lst_peptideDNASequence = []
lst_annotatedPeptideDNASequence = [["Gene Name",'Protein Name','Organism', 'Accession Number', 'Sequence of Interest']]
KSF_count = 0
ESF_count = 0
for entry in lst_s2J_results:

    uppercase_aaSequence = entry[1].upper()
    # print(uppercase_aaSequence)
    growingDNASequence = []

    for i in uppercase_aaSequence:
        # print(codonDictionary[i])
        growingDNASequence.append(codonDictionary[i])
        # print(growingDNASequence)
    joinedDNASequence = ''.join(growingDNASequence)

    ### This is where we put the restriction site removal lines
    matchKSF = re.search(restrictionKSF, joinedDNASequence)
    matchESF = re.search(restrictionESF, joinedDNASequence)
    if matchKSF:
        joinedDNASequence = joinedDNASequence.replace(restrictionKSF, fixedKSF)
        # print('KSF replaced')
        KSF_count += 1
    if matchESF:
        joinedDNASequence = joinedDNASequence.replace(restrictionESF, fixedESF)
        # print('ESF replaced')
        ESF_count += 1

    lst_peptideDNASequence.append(joinedDNASequence)
    lst_annotatedPeptideDNASequence.append([entry[0],entry[3],entry[4],entry[5],joinedDNASequence])
print('*** Replacing HindIII sites with next optimal codons ***')
print('KSF count =',KSF_count , '; ESF count =',ESF_count)
KSF_ESF_total_count = KSF_count + ESF_count
print('Total KSG and ESF count =',KSF_ESF_total_count)

print('*** Peptide to DNA - Codon Optimized for E. Coli and J to TAG ***')
print(lst_peptideDNASequence[:3])
print('# of peptide sequences converted to DNA:', len(lst_peptideDNASequence))
print('')
print(lst_annotatedPeptideDNASequence[:20])

### Added on Oct 21 ###
### Checking if Restriction Sites are present in peptide oligo sequences
site01_counter = 0
site02_counter = 0
for entry in lst_peptideDNASequence:
    # print(entry)
    match01 = re.search(restrictionSite01, entry)
    # print(match01)
    match02 = re.search(restrictionSite02, entry)
    if match01:
        site01_counter += 1
        # print(site01_counter)
    if match02:
        site02_counter += 1
        print(site02_counter)

        print(re.search(restrictionSite02, entry))
        print(entry)
        print(entry)
print('--- Checking if chosen Restriction Sites exist in peptide oligo sequences ---')
print('Restriction site 01 count=', site01_counter, '; Restriction site 02 count =', site02_counter)
print('')





lst_readyToGoOligoSequences = []
for peptideDNASequence in lst_peptideDNASequence:
    temp_oligo_lst = []
    temp_oligo_lst.append(startingNucleotides)
    temp_oligo_lst.append(peptideDNASequence)
    temp_oligo_lst.append(endingNucleotides)
    joinedOligoSequence= ''.join(temp_oligo_lst)
    lst_readyToGoOligoSequences.append(joinedOligoSequence)
print('Flanking Sequences added:')
print(lst_readyToGoOligoSequences[:3])
print('# of peptide sequences with flanking sequences:', len(lst_peptideDNASequence))

print('')
print(' --- ', inputPattern, ' ---')
print('total oligos before duplicates removed:',len(lst_readyToGoOligoSequences))
###Removing Duplicates
lst_readyToGoOligoSequencesDupRemoved = list(set(lst_readyToGoOligoSequences))
print('total oligos after duplicates removed:' ,len(lst_readyToGoOligoSequencesDupRemoved))
print('Duplicates removed:', len(lst_readyToGoOligoSequences) - len(lst_readyToGoOligoSequencesDupRemoved))

patternFileName = inputPattern.replace('.','x')

oligosFileName = time_str + '-' + patternFileName+ str(inputPlusMinus) + '-FormatedOligos.csv'
oligosFileNameAnnotated = time_str + '-' + patternFileName+ str(inputPlusMinus) + '-FormatedOligosAnnotated.csv'

# print(lst_readyToGoOligoSequencesDupRemoved)


##### File making Lines below
## Comment it out or in depending on needs.
with open(oligosFileName, 'w', newline='') as temp_Oligo:
    writer = csv.writer(temp_Oligo, delimiter=',')
    for row in lst_readyToGoOligoSequencesDupRemoved:
        writer.writerow([row])
shutil.move(oligosFileName, "generatedFiles/")

with open(oligosFileNameAnnotated, 'w', newline='') as temp_Oligo_Annotated:
    writer = csv.writer(temp_Oligo_Annotated)
    writer.writerows(lst_annotatedPeptideDNASequence)
shutil.move(oligosFileNameAnnotated, "generatedFiles/")




