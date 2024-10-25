import re
import csv
import time
import shutil
import numpy as np
import math




###dictionary for codons


codonDictionary = {
    'A':'GCG','R':'CGC','N':'AAC','D': 'GAT','C': 'TGC',
    'E': 'GAA','Q': 'CAG','G': 'GGC','H': 'CAT','I':'ATT',
    'L':'CTG','K': 'AAA','M': 'ATG','F': 'TTT','P': 'CCG',
    'S': 'AGC','T': 'ACC','W': 'TGG','Y': 'TAT','V': 'GTG',
    's': 'AGC','J': 'TAG'
}

print(codonDictionary.keys())

print(codonDictionary.values())


print(codonDictionary.get('A'))

aaSequenceSample = 'ttPsEsPRAQAtJRLJtAJcPtPKVQsRCss'

uppercase_aaSequence = aaSequenceSample.upper()
print(uppercase_aaSequence)
growingDNASequence = []

for i in uppercase_aaSequence:
    print(codonDictionary[i])
    growingDNASequence.append(codonDictionary[i])
    print(growingDNASequence)
joinedDNASequence = ''.join(growingDNASequence)
print(joinedDNASequence)