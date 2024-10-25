import re
import csv
import time
import shutil

import numpy as np
import math

reformatted_FASTA_file_name = "generatedFiles/CN reformatted PTM FASTA.csv"

with open(reformatted_FASTA_file_name) as FASTA_file:
    FASTA_data = FASTA_file.read()
# print(FASTA_data[:10000])
delimited_by_newLine = FASTA_data.split('\n')
reformatted_FASTA_List = []
for entry in delimited_by_newLine:
    entry.split(',')
    reformatted_FASTA_List.append(entry.split(','))
# print(delimited_by_newLine[:5])
print(reformatted_FASTA_List[:10])