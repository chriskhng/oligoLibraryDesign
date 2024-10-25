## - This contains the functions used for asSeq2Int.py for clarity
import re
import csv
import time
import shutil
import numpy as np
import math

from oligoLibraryDesignFunctions import FASTA_reformat

time_str = time.strftime("%Y%m%d-%H%M%S") ## - Used to prefix output files.

fakeFASTA = ['FAKE01', 'MMEEIDRFQDPAAASISDRDCDAREEKQRELARKGsTTxDDssLKNGSMGsPVNQQPKKNNVMARTRLVVPNKGYSsLDQsPDEKPLVALDtDsDDDFDMSRYSSSGYSSAEQINQDLNIQLLKDGYRLDEIPDDEDLDLIPPKSVNPTCMCCQATSSTACHI']#, ['Fake2', 'MNPSAAVIFCLILLGLSGTQGIPLARTVRCNCIHIDDGPVRMRAIGKLEIIPASLSCPRVEIIATMKKNDEQRCLNPESKTIKNLMKAFSQKRSKRA']]
user_input_SoI = "s..s..s"
lst_masterFASTA_SAMPLE = [['FAKE1', 'xxxxxABCDE123456789012s45s78s012345678901ABCDEyyyyy'], ['Cox7a2', 'MLRNLLALRQIAQRTISTTSRRHFENKVPEKQkLFQEDNGMPVHLkGGASDALLYRATMALTLGGTAYAIYLLAMAAFPKKQ'], ['FAKE02', 'xxxxxafdguauhfjnasdfs45s78swheifgsdgfbakjsdfhauisgiuafghluaig'], ['FAKE03', '123123ASDHIASHIUDsXXsYYsZZsGHBHUDSAFH123123'],  ['FAKE04', 'AAAAAAAAAAAAAAAAsAAsAAsBBBBBBCCCCCCCCsCCsCCsDDDDDDDDDDDDDDDDDDDD']]
lst_SoI_results = [['Gene Name', 'Sequence Hit', 'multiple?']]
pattern = 's..s..s'
inputPlusMinus = 15
spanPlusMinus = inputPlusMinus - math.ceil(len(pattern)/2)+1


FASTA = FASTA_reformat("Phosphosite_PTM_seq.fasta")

counter = 0
print('FOR LOOP BELOW')
for protein_entry in FASTA:
    # print(protein_entry)
    searchSpan = re.search(pattern, protein_entry[1])
    # print(result)
    if searchSpan:
        lowRange = np.clip(searchSpan.start() - (spanPlusMinus), 0, None)
        highRange = searchSpan.end() + (spanPlusMinus)
        print(protein_entry[1][lowRange:highRange])  ##THIS IS THE ONE
        lst_SoI_results.append([protein_entry[0],protein_entry[1][
              lowRange:highRange]])
        matches = re.findall(pattern, protein_entry[1])
        if len(matches) > 1:
            print('more than one')
            print(matches[1:])
            counter = counter + 1
    else:
        continue


print(lst_SoI_results[0:50])
print(len(lst_SoI_results))
print(counter)
# print("--- Searching for SoI in FASTA data ---")
# lst_SoI_results = [['Protein Name', 'Sequence Hit', 'Position of Hit']]




#
# ###BELOW IS IFROM ASSEQ2INT
# start_time = time.time()
# print("--- Searching for SoI in FASTA data ---")
# lst_SoI_results = [['Protein Name', 'Sequence Hit', 'Position of Hit']]
# ## - the user input SoI will be .upper, in case user input in lowercase
# upper_SoI = user_input_SoI.upper()
# ## - Then, if the user used x to denote "any", it will be replaced by .
# regex_SoI = upper_SoI.replace("X", ".")
# ## - the reformated user-input is then placed into RegEx pattern match syntax
# pattern = '^' + regex_SoI + "$"
# # SoI_check = regex_SoI.find("M")
# #
# # print(SoI_check)
#
# ## - The protein sequence in each protein entry in lst_masterFASTA will be iterated thru searching for SoI
# for protein_entry in lst_masterFASTA:
# ## protein_entry looks like:
#     #protein_entry = ['YAL001C', 'TFC3', 'SGDID:S000000001', 'MVLTIYPDELVQIVSDKI..."
#     protein_seq = protein_entry[3]
# ## - The following cuts each protein sequence into blocks with the length of the input SoI
#     for i in range(len(protein_seq)-(len(user_input_SoI)-1)):
#         protein_sequence_block =protein_seq[i:(i + len(user_input_SoI))]
#         # The following commented code checks to see if I cut the blocks correctly
#         # print(protein_seq[i:(i+len(pattern)-2)])
#         result = re.match(pattern, protein_sequence_block)
# ## - If a protein sequence block matches the user-input sequence, it will print:
#         if result: ############## - the print commands below may be omitted for the final code. keep for testing.
#             # print(upper_SoI + " found in:")
#             # print(protein_entry[1]) # the protein name and ...
#             # print("as " + protein_seq[i:(i + len(pattern) - 2)] + ' which starts at position: ' + str(i+1) + "\n")
#             # The following appends to a list that will be used to populate the output .txt file
#             lst_SoI_results.append([protein_entry[1], protein_seq[i:(i + len(pattern) - 2)], str(i+1)])
# ## - the actual sequence it matched, and which amino acid position the match started in the protein sequence
#         else:
#             continue
#
# ## - Here we define the output SeqHit file name
# SeqHitFileName = time_str +'SeqHit.csv'
# ## - lst_results is then written into SeqHitFileName with csv module
# with open(SeqHitFileName, 'w', newline='') as temp_SeqHit_file:
#     writer = csv.writer(temp_SeqHit_file)
#     writer.writerows(lst_SoI_results)
# shutil.move(SeqHitFileName, "SeqHitResults/")
#
# print("--- SoI search in FASTA took %s seconds ---" % (time.time() - start_time))
# print(str(len(lst_SoI_results) - 1)+ " number of SeqHits (this count includes "
#                                      "multiple hits within one protein sequence).")
# # return lst_SoI_results
# ## - #### The .py code is over for the purposes of the CIBR final project
# ## - #### Next, temp_results.csv will be copied and renamed to the date-time.csv
# ## - #### Finally, this renamed csv will be rsync'ed
# ## - #### (currently set to a temp_dir in the current dir, but can be changed to a local dir).