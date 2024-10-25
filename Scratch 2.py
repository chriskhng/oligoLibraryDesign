## - This contains the functions used for asSeq2Int.py for clarity
import re
import math
import numpy as np
import csv
import time
import shutil
#
# FASTA_file_name = "Phosphosite_PTM_seq.fasta"
#
# lst_masterFASTA = []
# arb_num_FASTA = 0
#
# while arb_num_FASTA == 0:
#     try:
#         start_time = time.time()
#         print("--- Reformatting FASTA entries ---")
#
#         with open(FASTA_file_name) as FASTA_file: #Actual code will ask for the name of the .FASTA
#             FASTA_data = FASTA_file.read()
#         ## - File can close now
#         FASTA_newlineM = FASTA_data.replace('\nM', '|M')
#         FASTA_newlineReplaced = FASTA_newlineM.replace('\n','')
#         FASTA_commaDelimited = FASTA_newlineReplaced.replace('|',',')
#         ## - lst_masterFASTAwithAuthors is a list of each protein entry, containing protein name, sequence, etc, in a single string
#         lst_masterFASTAwithAuthors = FASTA_commaDelimited.split('>GN:')   ## I split the starting FASTA with \n>GN as that's what's separtating the each entries in the FASTA:
#
#         ## - lst_rawFASTA is a list of each protein entry, containing protein name, sequence, etc, in a single string
#         lst_rawFASTA = lst_masterFASTAwithAuthors[1:-1] ##This removes the author watermark at the top which is index 0 here.
#         print(lst_masterFASTAwithAuthors[0:3])
#         print("/n ******************************************************")
#         print(lst_rawFASTA[0:3])
#         ## - relevant info will be pulled out from each protein entry
#         ## - We are interested in the protein name and protein sequence
#         for protein_entry in lst_rawFASTA:
#         ## - The protein name is pulled from each FASTA entry and saved into FASTAproteinNames
#         ## - This is done by indexing from 0 to the first comma (which delimits name from the next entry)
#             FASTAgeneNames = protein_entry[0: protein_entry.find(',')]
#             proteinSeqFinal = protein_entry[protein_entry.find(',M')+1: -1]
#
#         ## - Now, we put protein name and protein sequence together into a list for each entry
#         ## - First, the protein names are separated into genomic name, common name, and ID name
#             lst_FASTAentry = FASTAgeneNames.split(' ')
#         ## - Next, the protein sequence is appended to end of the names.
#             lst_FASTAentry.append (proteinSeqFinal)
#         ## - Finally, the list for the entry is appended to the growing list of protein entries from the FASTA.
#             lst_masterFASTA.append (lst_FASTAentry)
#
#
#         ## - At this point, lst_masterFASTA contains artifacts from the reformatting above.
#         ## - Namely, ">" and "*" remains on the first of first and last of last indices respectively.
#         ## - These artifacts are removed in the following 2 lines.
#         (lst_masterFASTA[0])[0] = (lst_masterFASTA[0])[0].replace('>', '')
#         (lst_masterFASTA[-1])[-1] = (lst_masterFASTA[-1])[-1].replace('*', '')
#         print("--- Reformatted FASTA entries in %s seconds ---" % (time.time() - start_time))
#         arb_num_FASTA = 1
#         print(lst_masterFASTA[0:3])
#
#     except FileNotFoundError:
#         FASTA_file_name = input ("The FASTA file name is invalid, try again: ")
#         ## - This list of sublists of protein entries will be used in Part 2.
#         ## - This list will be searched for the user-input SoI to output the name of proteins containing the SoI

# INITIAL PARAMATERS
lst_SoI_results = [['Gene Name', 'Sequence Hit']]
pattern = 's..s..s'
inputPlusMinus = 15
fakeFASTA = ['FAKE01', 'MMEEIDRFQDPAAASISDRDCDAREEKQRELARKGsTTxDDssLKNGSMGsPVNQQPKKNNVMARTRLVVPNKGYSsLDQsPDEKPLVALDtDsDDDFDMSRYSSSGYSSAEQINQDLNIQLLKDGYRLDEIPDDEDLDLIPPKSVNPTCMCCQATSSTACHI']#, ['Fake2', 'MNPSAAVIFCLILLGLSGTQGIPLARTVRCNCIHIDDGPVRMRAIGKLEIIPASLSCPRVEIIATMKKNDEQRCLNPESKTIKNLMKAFSQKRSKRA']]
# print(re.search('s..s..s', fakeFASTA[1]))



# string1 = "MMEEIDRFQDPAAASISDRDCDAREEKQRELARKGsTTsDDssLKNGSMGsPVNQQPKKNNVMARTRLVVPNKGYSsLDQsPDEKPLVALDtDsDDDFsDFsERsDMSRYSSSGYSSAEQINQDLNIQLLKDGYRLDEIPDDEDLDLIPPKSVNPTCMCCQATSSTACHI"
string1 = 'x56789012s45s78s012345678901ABCDEyyyyy'

print(re.search('s..s..s', string1))
print(re.findall('s..s..s', string1))
print(re.match('s..s..s', string1))
searchSpan = (re.search(pattern, string1))
# print(searchSpan)
print(searchSpan.start())
print(searchSpan.end())
# print('span:',searchSpan.span())
# print('start -12:',searchSpan.start()-12)
# print('end +12:',searchSpan.end()+12)
# print(string1[0:])
print('XXX' ,string1[searchSpan.start()-12:searchSpan.end()+12]) ##THIS IS THE ONE
# print('pattern length:', len(pattern))
# print('pattern length/2:', (len(pattern)/2))
# print('pattern length/2 ROUND:', math.ceil(len(pattern)/2))
# print('15 - pattern length/2 ROUND:', inputPlusMinus - math.ceil(len(pattern)/2))
print('15 - pattern length/2 ROUND +1 :', inputPlusMinus - math.ceil(len(pattern)/2)+1) ###this is the one too

spanPlusMinus = inputPlusMinus - math.ceil(len(pattern)/2)+1

print('spanPlusMinus:', spanPlusMinus)

# print(string1[searchSpan.start()-(inputPlusMinus - math.ceil(len(pattern)/2)+1):searchSpan.end()+(inputPlusMinus - math.ceil(len(pattern)/2)+1)])
print('the correct answer before for loop efforts:',string1[searchSpan.start()-(spanPlusMinus):searchSpan.end()+(spanPlusMinus)]) ##THIS IS THE ONE
print( searchSpan.start())
print( searchSpan.start()-(spanPlusMinus))
print((np.clip(searchSpan.start()-(spanPlusMinus),0,None)))
clipTest = searchSpan.start()-(spanPlusMinus)
print(clipTest)
print(np.clip(clipTest,0,10))

lst_masterFASTA_SAMPLE = [['FAKE1', 'xxxxxABCDE123456789012s45s78s012345678901ABCDEyyyyy'], ['Cox7a2', 'MLRNLLALRQIAQRTISTTSRRHFENKVPEKQkLFQEDNGMPVHLkGGASDALLYRATMALTLGGTAYAIYLLAMAAFPKKQ'], ['FAKE02', 'xxxxxafdguauhfjnasdfs45s78swheifgsdgfbakjsdfhauisgiuafghluaig'], ['FAKE03', '123123ASDHIASHIUDsXXsYYsZZsGHBHUDSAFH123123'],  ['FAKE04', 'AAAAAAAAAAAAAAAAsAAsAAsBBBBBBCCCCCCCCsCCsCCsDDDDDDDDDDDDDDDDDDDD']]
#xxxxxABCDE123456789012s45s78s012345678901ABCDEyyyyy
# print(lst_masterFASTA_SAMPLE)
# print(lst_masterFASTA_SAMPLE[0])
# print(lst_masterFASTA_SAMPLE[0][0])
# print(lst_masterFASTA_SAMPLE[0][1])

print(lst_masterFASTA_SAMPLE[0][1][searchSpan.start()-(spanPlusMinus):searchSpan.end()+(spanPlusMinus)]) ##THIS IS THE ONE
print(lst_masterFASTA_SAMPLE[0][1][np.clip(searchSpan.start()-(spanPlusMinus),0,None):searchSpan.end()+(spanPlusMinus)]) ##THIS IS THE ONE
print(searchSpan.start()-(spanPlusMinus))
# print(lst_masterFASTA_SAMPLE[1][1])
# print(lst_masterFASTA_SAMPLE[1][1])



print('FOR LOOP BELOW')
for protein_entry in lst_masterFASTA_SAMPLE:
    # print(protein_entry)
    searchSpan = re.search(pattern, protein_entry[1])
    # print(result)
    if searchSpan:
        lowRange = np.clip(searchSpan.start() - (spanPlusMinus), 0, None)
        highRange = searchSpan.end() + (spanPlusMinus)
        print(protein_entry[1][lowRange:highRange])  ##THIS IS THE ONE
        lst_SoI_results.append([protein_entry[0],protein_entry[1][
              lowRange:highRange]])
    else:
        continue






print(lst_SoI_results)


string1 = 'x56789012s45s78s012345678901ABCDEyyyyy'
string2 = ''
FakeFASTA = ['FAKE01', 'MMEEIDRFQDPAAASISDRDCDAREEKQsREsLAsRKsGsTTxDDssLKNGSMGsPVNQQPKKNsNVsMAsRTRLVVPNKGLDLIPPKSVNPTCMCCQATSSTACHI'], ['Fake2', 'MNPSAAVIFCLILLGLSGTQGIPLARTVRCNCIHIDDGPVRMRAIGKLEIIPASLSCPRVEIIATMKKNDEQRCLNPESKTIKNLMKAFSQKRSKRA']
print(FakeFASTA[0][1])
print('research Fake01:' , re.search('s..s..s', FakeFASTA[0][1]))
print('re.findall Fake01:' , re.findall('s..s..s', FakeFASTA[0][1]))
print('re.findall Fake01:' , re.findall('s..s..s', FakeFASTA[0][1]))
matches = re.findall('s..s..s', FakeFASTA[0][1])
print('matches:', matches)
print('matches 1:', matches[1])
print(re.match('s..s..s', string1))
print(len(matches))
# counter = 0
# if len(matches) > 1:
#     print('more than one')
#     counter = counter + 1
# else:
#     print('just one')
# print(counter)
# # # #
# # #
# # #
# # #
# # a = -1
# #
# # print(a)
# # print(np.clip(a,2,8))
#
lst_SoI_resultsB = [['Gene Name', 'Sequence Hit']]
FindallFASTATEST = [['FAKE01', '12sddsxxs901asdasdasd23456sAAsBBs45678901234567ssssssssssssss231231ssGsLGs301'], ['Fake2', 'xxxxxABCDE123456789012s45s78s012345678901ABCDEyyyyy']]
matchesTest =  ['sddsxxs','sAAsBBs', 'ssGsLGs']

# print(re.search(matchesTest[0],FindallFASTATEST[0][1]))

print('FOR LOOP BELOW')
for protein_entry in FindallFASTATEST:
    # print(protein_entry)
    searchSpan = re.search(pattern, protein_entry[1])
    # print(searchSpan)
    if searchSpan:
        matches = re.findall(pattern, protein_entry[1])
        if len(matches) > 1:
            # print('more than one')
            print(matches)
            for multiPattern in matches:
                multiSearchSpan = re.search(multiPattern, protein_entry[1])
                # print(multiSearchSpan)
                print(multiPattern)
                multiLowRange = np.clip(multiSearchSpan.start() - (spanPlusMinus), 0, None)
                multiHighRange = multiSearchSpan.end() + (spanPlusMinus)
                lst_SoI_resultsB.append([protein_entry[0], protein_entry[1][
                                                          multiLowRange:multiHighRange]])
        else:
            lowRange = np.clip(searchSpan.start() - (spanPlusMinus), 0, None)
            highRange = searchSpan.end() + (spanPlusMinus)
            # print(protein_entry[1][lowRange:highRange])  ##THIS IS THE ONE
            lst_SoI_resultsB.append([protein_entry[0],protein_entry[1][
                                                      lowRange:highRange]])

    else:
        continue
print(lst_SoI_resultsB)


print(re.findall('s.s.s', 's1s2s3s4s5s6s7s8s9s0ssssssssss'))
print(re.findall('s.s.s', 'B1s2s3s4s5s6s7s8s9s0ssssssssss'))