import re
import csv
import time
import shutil
import numpy as np
import math

from oligoLibraryDesign import joinedDNASequence

### index 0 is a negative. index 1 is a positive for 01
### Index 3 is a positiv efor 02. Index 4 is a positive for both.
### Index 5 is positive for 01
lst_peptideDNASequence = ['AGCCGCCTGCCGAAATGGAGCGAAGGCCCGTGCGGCTAGCCGGGCTAGCCGGGCTAGAACATTGAAGGCTTTGTGTGCAAATTTAGCTTTAAA', 'ACCACCCCGAGCGAAAGCCCGCGCGCGCAGGCGAACCAAGCCTAGCGCCTGTAGACCGCGTAGTGCCCGACCCCGAAAGTGCAGAGCCGCTGCAGCAGC', 'ATGGTGGCGGCGAAAAAAACCAAAAAATAGCTGGAATAGATTAACTAGCGCCTGCAGCTGGTGAAAGCTTTGAAAAGCGGCAAATATGTG', 'AGCCGCCTGCCGAAACCAAGATGGAGCGAAGGCCCGTGCGGCTAGCCGGGCTAGCCGGGCTAGAACATTGAAGGCTTTGTGAAGCTTTGCAAATTTAGCTTTAAA', 'ACCACCCCGAGCGAAAGCCCGCGCGCGCAGGCGAACCAAGCCTAGCGCCTGTAGACCGCGTAGTGCCCGACCCCGAAAGTGCAGAGCCGCTGCAGCAGC']

restrictionSite01 = 'ACCAAG' ###DpnI
restrictionSite02 = 'AAGCTT' ###HindIII


print('findall print readout for Site01')
print(re.findall(restrictionSite01, lst_peptideDNASequence[0]))
print(re.findall(restrictionSite01, lst_peptideDNASequence[1]))
print(re.findall(restrictionSite01, lst_peptideDNASequence[2]))
print(re.findall(restrictionSite01, lst_peptideDNASequence[3]))
print('')

site01_counter = 0
site02_counter = 0
for entry in lst_peptideDNASequence:
    print(entry)
    match01 = re.search(restrictionSite01, entry)
    # print(match01)
    match02 = re.search(restrictionSite02, entry)

    if match01:
        site01_counter += 1
        # print(site01_counter)
    if match02:
        site02_counter += 1
        print(re.search(restrictionSite02, entry))
        # print(site02_counter)
print('--- Checking if chosen Restriction Sites exist in peptide oligo sequences ---')
print('Restriction site 01 count=', site01_counter)
print('Restriction site 02 count =', site02_counter)



print('')

print('findall print readout for Site02')
print(re.findall(restrictionSite02, lst_peptideDNASequence[0]))
print(re.findall(restrictionSite02, lst_peptideDNASequence[1]))
print(re.findall(restrictionSite02, lst_peptideDNASequence[2]))
print(re.findall(restrictionSite02, lst_peptideDNASequence[3]))

### index 0 is a negative. index 1 is a positive for KSF
### Index 2 is a positive for ESF.

lst_peptideDNASequence = ['AGCCGCCTGCCGAAATGGAGCGAAGGCCCGTGCGGCTAGCCGGGCTAGCCGGGCTAGAACATTGAAGGCTTTGTGTGCAAATTTAGCTTTAAA', 'CCGACCTTTGGCAAAAGCTTTCATTTTGATCCGCTGTAGAGCGGCTAGCGCAGCTAGAGCCTGAAAAGCGCGCAGGGCACCGGCTTTGAACTG', 'GCGGAAGCGGATCATAGCGGCGGCAGCGATCGCAACTAGATGGATTAGGTGGATTAGTGCTGCAGCCTGAAAAAAACCGAAAGCTTTCAGAAC']


print('')
print('*******************************')
print('*** Find and replace lines: ***')
print('*******************************')
print('')

site01_counter = 0
site02_counter = 0
for entry in lst_peptideDNASequence:
    print(entry)
    match01 = re.search(restrictionSite01, entry)
    # print(match01)
    match02 = re.search(restrictionSite02, entry)

    if match01:
        site01_counter += 1
        # print(site01_counter)
    if match02:
        site02_counter += 1
        print(re.search(restrictionSite02, entry))
        # print(site02_counter)
print('--- Checking if chosen Restriction Sites exist in peptide oligo sequences ---')
print('Restriction site 01 count=', site01_counter)
print('Restriction site 02 count =', site02_counter)
print('')

joinedDNASequence = 'CCGACCTTTGGCAAAAGCTTTCATTTTGATCCGCTGTAGAGCGGCTAGCGCAGCTAGAGCCTGAAAAGCGCGCAGGGCACCGGCTTTGAACTG'  ##with KSF
joinedDNASequence = 'GGCGGCGGCTTTTTTGAAAGCTTTAAACGCGTGATTCGCTAGCGCTAGCAGTAGATGGATGCGATGGGCCTGAGCAACAAAAAACCGAACACC' ## with ESF
restrictionKSF = 'AAAAGCTTT'
restrictionESF = 'GAAAGCTTT'

fixedKSF = 'AAAAGTTTT'
fixedESF = 'GAAAGTTTT'

## Solution 1

KSF_count = 0
ESF_count = 0
print('solution 1')

print('original sequence =')
print(joinedDNASequence)
matchKSF = re.search(restrictionKSF,joinedDNASequence)

matchESF = re.search(restrictionESF,joinedDNASequence)
if matchKSF:
    joinedDNASequence = joinedDNASequence.replace(restrictionKSF, fixedKSF)
    # print('KSF replaced')
    KSF_count += 1
if matchESF:
    joinedDNASequence = joinedDNASequence.replace(restrictionESF, fixedESF)
    # print('ESF replaced')
    ESF_count += 1


print('KSF count =',KSF_count)
print('ESF count =',ESF_count)
KSF_ESF_total_count = KSF_count + ESF_count
print('Total KSG and ESF count =',KSF_ESF_total_count)

print('new:')
print(joinedDNASequence)