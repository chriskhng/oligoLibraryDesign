import re
# import csv
# import time
# import shutil
# import numpy as np
# import math
### Inputs ###
pattern = 's..s..s'
inputPlusMinus = 20

###standard length####
###standard length####
standardLength = inputPlusMinus*2 + len(pattern)%2 #%2 will give 1 when pattern length is odd. 0 when pattern length is even.
# print(standardLength)

####Difference from standard length



SampleSoIOutput = [['CD93', 'QPLLPsRLPKWSEGPCGsPGsPGsNIEGFVCKFsFKGMCRP', 'sPGsPGs'], ['GPHN', 'TASLSttPsEsPRAQAtsRLstAscPtPKVQsRCsskENIL', 'sRLstAs'], ['Rpl30', 'MVAAKKTKksLEsINsRLQLVMkSGkYVLGYkQ', 'sLEsINs'], ['Pdap1', 'QEEGGDGAsGDPkkEKKsLDsDEsEDEDDDyQQKRKGVEGL', 'sLDsDEs'], ['Cmklr1', 'ALFSRLANALsEDtGPssyPsHRsFtKMssLNEkAsVNEKE', 'syPsHRs'], ['Arpp21', 'EREEEYQRVRERIFAHDsVCsQEsLFLDNSRLQEDMHICNE', 'sVCsQEs'], ['Neb', 'LPQQRsSsVAtQQTTLSsIPsHPstAGKIFRAIYDYIAADA', 'sIPsHPs'], ['Ryr2', 'RTREGDsMALYNRTRRIsQtsQVsIDAAHGYsPRAIDMSNV', 'sQtsQVs'], ['Sub1', 'MPKSKELVsssssGsDsDsEVEKKLKRKKQAV', 'sssssGs'], ['Meaf6', 'SkNDRRNRKFKEAERLFsKssVTsAAAVSALAGVQDQLIEK', 'sKssVTs'], ['MACF1', 'PSsRAAsPTRsssSAsQsNHsCTsMPssPAtPASGTKVIPS', 'sNHsCTs'], ['PCLO', 'STVQLAPsPPKsPKVLYsPIsPLsPGKALESAFVPYEKPLP', 'sPIsPLs'], ['BCL11A', 'GEsASGGLskkLLLGsPssLsPFskRIKLEKEFDLPPAAMP', 'ssLsPFs'], ['Pgls', 'RAASCLEGDRGRFALGLsGGsLVsMLARDLPAAAAPAGPAS', 'sGGsLVs'], ['Cacng2', 'AsAITRIPsyRYRYQRRsRsssRsTEPsHsRDAsPVGVKGF', 'sRsssRs'], ['ANK3', 'sLPAEGyMGFsLGARSAsLRsFssDrsYTLNRSSYARDsMM', 'sLRsFss'], ['MGP', 'LLLSILAALAVAALCYEsHEsLEsYEINPFINRRNANSFIS', 'sHEsLEs'], ['KMT2A', 'KPRGRPRsGsDRNsAILsDPsVFsPLNKsEtKsGDKIKKKD', 'sDPsVFs'], ['KMT2A', 'QsPENESNDRRSRRysVsERsFGsRTTKKLSTLQSAPQQQT', 'sERsFGs']]

singleSoISampleEntry = [['DoubleFlank', 'QPLLPsRLPKWSEGPCGsPGsPGsNIEGFVCKFsFKGMCRP']]


###THIS IS IT
patternSoI = SampleSoIOutput[1][2]
extendedSequence = SampleSoIOutput[1][1]

for entry in SampleSoIOutput:
    print("entry:" ,entry)
    print(re.search(entry[2],entry[1])) ## This searchs for the sequence hit in the plus mins sequence
    firstSerineIndex = (re.search(entry[2],entry[1]).start())  ### first serine
    print('firstSerineIndex:',firstSerineIndex)
    thirdSerineIndex=(re.search(entry[2],entry[1]).end()) -1 #### Third serine
    print('thirdSerineIndex:',thirdSerineIndex)
    ###central serine
    # print('re.research output with [1:]:', re.search('s', pattern[1:])) # re.findall output: <re.Match object; span=(2, 3), match='s'>
    centralSerineIndexTemp = re.search('s', pattern[1:])
    # print('centralSerineIndexTemp:', centralSerineIndexTemp)
    # print('centralSerineIndexTemp.start():', centralSerineIndexTemp.start())
    centralSerineIndexAdd = centralSerineIndexTemp.start() + 1
    # print('centralSerineAdd:', centralSerineIndexAdd)
    centralSerineIndex = centralSerineIndexAdd +firstSerineIndex
    print('centralSerineIndex:', centralSerineIndex)

    seqIntoList = list(entry[1])
    print(seqIntoList)
    seqIntoList[firstSerineIndex] = ("J")
    seqIntoList[centralSerineIndex] = ("J")
    seqIntoList[thirdSerineIndex] = ("J")
    print(seqIntoList)
    # print(tempList[17])
    rejoinedList = ''.join(seqIntoList)
    print(rejoinedList)
    print('')