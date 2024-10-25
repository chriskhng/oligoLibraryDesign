import re
# import csv
# import time
# import shutil
# import numpy as np
# import math
### Inputs ###
pattern = 's..s..s'
inputPlusMinus = 20

#####This scratch is to figure out how swap the target s with J


###standard length####
###standard length####
standardLength = inputPlusMinus*2 + len(pattern)%2 #%2 will give 1 when pattern length is odd. 0 when pattern length is even.
print(standardLength)

####Difference from standard length





# SampleSoIOutput = [['Gene Name', 'Sequence Hit'], ['DoubleFlank', 'QPLLPsRLPKWSEGPCGsPGsPGsNIEGFVCKFsFKGMCRP'], ['GPHN', 'TASLSttPsEsPRAQAtsRLstAscPtPKVQsRCsskENIL'], ['Shortstart8', 'MVAAKKTKksLEsINsRLQLVMkSGkYVLGYkQ'], ['Pdap1', 'QEEGGDGAsGDPkkEKKsLDsDEsEDEDDDyQQKRKGVEGL'], ['Cmklr1', 'ALFSRLANALsEDtGPssyPsHRsFtKMssLNEkAsVNEKE'], ['Arpp21', 'EREEEYQRVRERIFAHDsVCsQEsLFLDNSRLQEDMHICNE'], ['Neb', 'LPQQRsSsVAtQQTTLSsIPsHPstAGKIFRAIYDYIAADA'], ['Ryr2', 'RTREGDsMALYNRTRRIsQtsQVsIDAAHGYsPRAIDMSNV'], ['Sub1', 'MPKSKELVsssssGsDsDsEVEKKLKRKKQAV'], ['Meaf6', 'SkNDRRNRKFKEAERLFsKssVTsAAAVSALAGVQDQLIEK'], ['MACF1', 'PSsRAAsPTRsssSAsQsNHsCTsMPssPAtPASGTKVIPS'], ['PCLO', 'STVQLAPsPPKsPKVLYsPIsPLsPGKALESAFVPYEKPLP'], ['BCL11A', 'GEsASGGLskkLLLGsPssLsPFskRIKLEKEFDLPPAAMP'], ['Pgls', 'RAASCLEGDRGRFALGLsGGsLVsMLARDLPAAAAPAGPAS'], ['Cacng2', 'AsAITRIPsyRYRYQRRsRsssRsTEPsHsRDAsPVGVKGF'], ['ANK3', 'sLPAEGyMGFsLGARSAsLRsFssDrsYTLNRSSYARDsMM'], ['MGP', 'LLLSILAALAVAALCYEsHEsLEsYEINPFINRRNANSFIS'], ['KMT2A', 'KPRGRPRsGsDRNsAILsDPsVFsPLNKsEtKsGDKIKKKD'], ['KMT2A', 'QsPENESNDRRSRRysVsERsFGsRTTKKLSTLQSAPQQQT']]
# SampleSoIOutput = [['Gene Name', 'Sequence Hit'], ['DoubleFlank', 'QPLLPsRLPKWSEGPCGsPGsPGsNIEGFVCKFsFKGMCRP', 'sPGsPGs'], ['GPHN', 'TASLSttPsEsPRAQAtsRLstAscPtPKVQsRCsskENIL', 'sRLstAs'], ['Shortstart8', 'MVAAKKTKksLEsINsRLQLVMkSGkYVLGYkQ', 'sLEsINs'], ['Pdap1', 'QEEGGDGAsGDPkkEKKsLDsDEsEDEDDDyQQKRKGVEGL', 'sLDsDEs'], ['Cmklr1', 'ALFSRLANALsEDtGPssyPsHRsFtKMssLNEkAsVNEKE', 'syPsHRs'], ['Arpp21', 'EREEEYQRVRERIFAHDsVCsQEsLFLDNSRLQEDMHICNE', 'sVCsQEs'], ['Neb', 'LPQQRsSsVAtQQTTLSsIPsHPstAGKIFRAIYDYIAADA', 'sIPsHPs'], ['Ryr2', 'RTREGDsMALYNRTRRIsQtsQVsIDAAHGYsPRAIDMSNV', 'sQtsQVs'], ['Sub1', 'MPKSKELVsssssGsDsDsEVEKKLKRKKQAV', 'sssssGs'], ['Meaf6', 'SkNDRRNRKFKEAERLFsKssVTsAAAVSALAGVQDQLIEK', 'sKssVTs'], ['MACF1', 'PSsRAAsPTRsssSAsQsNHsCTsMPssPAtPASGTKVIPS', 'sNHsCTs'], ['PCLO', 'STVQLAPsPPKsPKVLYsPIsPLsPGKALESAFVPYEKPLP', 'sPIsPLs'], ['BCL11A', 'GEsASGGLskkLLLGsPssLsPFskRIKLEKEFDLPPAAMP', 'ssLsPFs'], ['Pgls', 'RAASCLEGDRGRFALGLsGGsLVsMLARDLPAAAAPAGPAS', 'sGGsLVs'], ['Cacng2', 'AsAITRIPsyRYRYQRRsRsssRsTEPsHsRDAsPVGVKGF', 'sRsssRs'], ['ANK3', 'sLPAEGyMGFsLGARSAsLRsFssDrsYTLNRSSYARDsMM', 'sLRsFss'], ['MGP', 'LLLSILAALAVAALCYEsHEsLEsYEINPFINRRNANSFIS', 'sHEsLEs'], ['KMT2A', 'KPRGRPRsGsDRNsAILsDPsVFsPLNKsEtKsGDKIKKKD', 'sDPsVFs'], ['KMT2A', 'QsPENESNDRRSRRysVsERsFGsRTTKKLSTLQSAPQQQT', 'sERsFGs']]
SampleSoIOutput = [['Gene Name', 'Extended Sequence', 'Sequence Hit'], ['CD93', 'QPLLPsRLPKWSEGPCGsPGsPGsNIEGFVCKFsFKGMCRP', 'sPGsPGs'], ['GPHN', 'TASLSttPsEsPRAQAtsRLstAscPtPKVQsRCsskENIL', 'sRLstAs'], ['Rpl30', 'MVAAKKTKksLEsINsRLQLVMkSGkYVLGYkQ', 'sLEsINs'], ['Pdap1', 'QEEGGDGAsGDPkkEKKsLDsDEsEDEDDDyQQKRKGVEGL', 'sLDsDEs'], ['Cmklr1', 'ALFSRLANALsEDtGPssyPsHRsFtKMssLNEkAsVNEKE', 'syPsHRs'], ['Arpp21', 'EREEEYQRVRERIFAHDsVCsQEsLFLDNSRLQEDMHICNE', 'sVCsQEs'], ['Neb', 'LPQQRsSsVAtQQTTLSsIPsHPstAGKIFRAIYDYIAADA', 'sIPsHPs'], ['Ryr2', 'RTREGDsMALYNRTRRIsQtsQVsIDAAHGYsPRAIDMSNV', 'sQtsQVs'], ['Sub1', 'MPKSKELVsssssGsDsDsEVEKKLKRKKQAV', 'sssssGs'], ['Meaf6', 'SkNDRRNRKFKEAERLFsKssVTsAAAVSALAGVQDQLIEK', 'sKssVTs'], ['MACF1', 'PSsRAAsPTRsssSAsQsNHsCTsMPssPAtPASGTKVIPS', 'sNHsCTs'], ['PCLO', 'STVQLAPsPPKsPKVLYsPIsPLsPGKALESAFVPYEKPLP', 'sPIsPLs'], ['BCL11A', 'GEsASGGLskkLLLGsPssLsPFskRIKLEKEFDLPPAAMP', 'ssLsPFs'], ['Pgls', 'RAASCLEGDRGRFALGLsGGsLVsMLARDLPAAAAPAGPAS', 'sGGsLVs'], ['Cacng2', 'AsAITRIPsyRYRYQRRsRsssRsTEPsHsRDAsPVGVKGF', 'sRsssRs'], ['ANK3', 'sLPAEGyMGFsLGARSAsLRsFssDrsYTLNRSSYARDsMM', 'sLRsFss'], ['MGP', 'LLLSILAALAVAALCYEsHEsLEsYEINPFINRRNANSFIS', 'sHEsLEs'], ['KMT2A', 'KPRGRPRsGsDRNsAILsDPsVFsPLNKsEtKsGDKIKKKD', 'sDPsVFs'], ['KMT2A', 'QsPENESNDRRSRRysVsERsFGsRTTKKLSTLQSAPQQQT', 'sERsFGs']]

singleSoISampleEntry = [['DoubleFlank', 'QPLLPsRLPKWSEGPCGsPGsPGsNIEGFVCKFsFKGMCRP']]


###THIS IS IT
patternSoI = SampleSoIOutput[1][2]
extendedSequence = SampleSoIOutput[1][1]

print(re.search(SampleSoIOutput[1][2],SampleSoIOutput[1][1])) ## This searchs for the sequence hit in the plus mins sequence
firstSerineIndex = (re.search(SampleSoIOutput[1][2],SampleSoIOutput[1][1]).start())  ### first serine
print('firstSerineIndex:',firstSerineIndex)
thirdSerineIndex=(re.search(SampleSoIOutput[1][2],SampleSoIOutput[1][1]).end()) -1 #### Third serine
print('thirdSerineIndex:',thirdSerineIndex)
###central serine
print('re.research output with [1:]:', re.search('s', pattern[1:])) # re.findall output: <re.Match object; span=(2, 3), match='s'>
centralSerineIndexTemp = re.search('s', pattern[1:])
print('centralSerineIndexTemp:', centralSerineIndexTemp)
print('centralSerineIndexTemp.start():', centralSerineIndexTemp.start())
centralSerineIndexAdd = centralSerineIndexTemp.start() + 1
print('centralSerineAdd:', centralSerineIndexAdd)
centralSerineIndex = centralSerineIndexAdd +firstSerineIndex
print('centralSerineIndex:', centralSerineIndex)

tempList = list(extendedSequence)
print(tempList)
tempList[firstSerineIndex] = ("J")
tempList[centralSerineIndex] = ("J")
tempList[thirdSerineIndex] = ("J")
print(tempList)
print(tempList[17])
rejoinedList = ''.join(tempList)
print(rejoinedList)

print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      ''
      ''
      'XXXXXXXXXXXXXXXXXXXXXXX')


###figuroing out how to swap target s to J
sampleSoI = singleSoISampleEntry[0][1] # Output: QPLLPsRLPKWSEGPCGsPGsPGsNIEGFVCKFsFKGMCRP
print(sampleSoI)

print(re.findall('s..s..s',sampleSoI))
matches = re.findall('s..s..s', sampleSoI)
print(matches) #output is ['sPGsPGs']
print('matches[0]:', matches[0]) #output is sPGsPGs
print(re.match('s..s..s', sampleSoI)) #Output is None
print(re.search('s..s..s', sampleSoI)) #Output is <re.Match object; span=(17, 24), match='sPGsPGs'>
searchResult = re.search('s..s..s', sampleSoI)
print(searchResult.start(),searchResult.end()) #outout is 17 24

print('sampleSoI:',sampleSoI)
print('sampleSoI[20]:', sampleSoI[20])




#firstSerineIndex = inputPlusMinus + len(pattern)%2
firstSerineIndex = searchResult.start()
print('firstSerineIndex:', firstSerineIndex) #17

###third Serine
thirdSerineIndex = searchResult.end()
print('thirdSerineIndex:', thirdSerineIndex)  #24

###central serine
print('re.research output with [1:]:', re.search('s', pattern[1:])) # re.findall output: <re.Match object; span=(2, 3), match='s'>
centralSerineIndexTemp = re.search('s', pattern[1:])
print('centralSerineIndexTemp:', centralSerineIndexTemp)
print('centralSerineIndexTemp.start():', centralSerineIndexTemp.start())
centralSerineIndex = centralSerineIndexTemp.start() + 1
print('centralSerineIndex:', centralSerineIndex)


tempList = list(sampleSoI)
print(tempList)
tempList[17] = ("J")
print(tempList)
print(tempList[17])
rejoinedList = ''.join(tempList)
print(rejoinedList)



