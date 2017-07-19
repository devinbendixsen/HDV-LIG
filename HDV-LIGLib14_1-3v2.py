# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 02:31:44 2016

@author: devinbendixsen
"""

import pickle
import numpy
mutnumHDVLIGLib14=pickle.load(open("mutnumHDVLIGLib14.p", "rb"))
binarysequences=pickle.load(open("binarysequencesHDVLIGLib14.p", "rb"))

#HDV---------------------------------------------------------------------------------
ribozyme_start=('GACTCCCA', 'GAATCCCA')
rep=1
while rep!=4:
    cleaved={}
    uncleaved={}
    PercentCleaved={}
    for sequence in open("HDV"+str(rep)+".fastq"):
        sequence=sequence.rstrip()
        if sequence.startswith('A') or sequence.startswith('C') or sequence.startswith('T') or sequence.startswith('G'):
            if 'TCCCATTAG' in sequence and 'GCGGCGGGAGTTG' in sequence and 'AGGGAGGAA' in sequence and 'CCGCCTCCT' in sequence and ('GACTCCCA' in sequence or 'GAATCCCA' in sequence) and 'CTATAGGA' not in sequence:#CTATAGGA is 3nt of the TATA box and 5nt Rz sequence
                if '7' not in sequence and ')' not in sequence and 'N' not in sequence and '<' not in sequence and 'AAAAAA' not in sequence:         
                    if 'GACCATTC' in sequence:
                        while not sequence.startswith(ribozyme_start):
                            sequence=sequence.replace(' ', '')[1:]
                        sequence=sequence[2:3]+sequence[12:13]+sequence[17:18]+sequence[27:28]+sequence[41:42]+sequence[44:45]+sequence[54:55]+sequence[57:58]+sequence[61:63]+sequence[67:68]+sequence[72:73]+sequence[74:76]       
                        if (sequence[:1]=='A' or sequence[:1]=='C') and (sequence[1:2]=='A' or sequence[1:2]=='G') and (sequence[2:3]=='G' or sequence[2:3]=='T') and (sequence[3:4]=='G' or sequence[3:4]=='C') and (sequence[4:5]=='G' or sequence[4:5]=='C') and (sequence[5:6]=='G' or sequence[5:6]=='T') and (sequence[6:7]=='G' or sequence[6:7]=='C') and (sequence[7:8]=='C' or sequence[7:8]=='T') and (sequence[8:9]=='C' or sequence[8:9]=='T') and (sequence[9:10]=='C' or sequence[9:10]=='T') and (sequence[10:11]=='G' or sequence[10:11]=='A') and (sequence[11:12]=='G' or sequence[11:12]=='C') and (sequence[12:13]=='A' or sequence[12:13]=='C') and (sequence[13:14]=='G' or sequence[13:14]=='C'):        
                            if sequence in uncleaved:
                                uncleaved[sequence]+=1
                            elif sequence not in uncleaved:
                                uncleaved[sequence]=1   
                    elif 'GACCATTC' not in sequence:
                        while not sequence.startswith(ribozyme_start):
                            sequence=sequence.replace(' ', '')[1:]
                        sequence=sequence[2:3]+sequence[12:13]+sequence[17:18]+sequence[27:28]+sequence[41:42]+sequence[44:45]+sequence[54:55]+sequence[57:58]+sequence[61:63]+sequence[67:68]+sequence[72:73]+sequence[74:76]       
                        if (sequence[:1]=='A' or sequence[:1]=='C') and (sequence[1:2]=='A' or sequence[1:2]=='G') and (sequence[2:3]=='G' or sequence[2:3]=='T') and (sequence[3:4]=='G' or sequence[3:4]=='C') and (sequence[4:5]=='G' or sequence[4:5]=='C') and (sequence[5:6]=='G' or sequence[5:6]=='T') and (sequence[6:7]=='G' or sequence[6:7]=='C') and (sequence[7:8]=='C' or sequence[7:8]=='T') and (sequence[8:9]=='C' or sequence[8:9]=='T') and (sequence[9:10]=='C' or sequence[9:10]=='T') and (sequence[10:11]=='G' or sequence[10:11]=='A') and (sequence[11:12]=='G' or sequence[11:12]=='C') and (sequence[12:13]=='A' or sequence[12:13]=='C') and (sequence[13:14]=='G' or sequence[13:14]=='C'):        
                            if sequence in cleaved:
                                cleaved[sequence]+=1
                            elif sequence not in cleaved:
                                cleaved[sequence]=1
    for key in cleaved:
        if key in uncleaved:
            PercentCleaved[key]=float(cleaved[key])/(float(uncleaved[key])+float(cleaved[key]))
        else:
            uncleaved[key]=0
            PercentCleaved[key]=float(cleaved[key])/(float(uncleaved[key])+float(cleaved[key]))
    print ("Number Cleaved "+str(rep), ":", len(cleaved))
    print ("Number Uncleaved "+str(rep), ":", len(uncleaved))
    print ("Number %Cleaved "+str(rep), ":", len(PercentCleaved))
    pickle.dump(PercentCleaved, open("PercentCleaved"+str(rep)+".p", "wb"))
    pickle.dump(cleaved, open("cleaved"+str(rep)+".p", "wb"))
    pickle.dump(uncleaved, open("uncleaved"+str(rep)+".p", "wb"))
    rep+=1 

PercentCleaved1=pickle.load(open("PercentCleaved1.p", "rb"))
PercentCleaved2=pickle.load(open("PercentCleaved2.p", "rb"))
PercentCleaved3=pickle.load(open("PercentCleaved3.p", "rb"))
cleaved1=pickle.load(open("cleaved1.p", "rb"))
cleaved2=pickle.load(open("cleaved2.p", "rb"))
cleaved3=pickle.load(open("cleaved3.p", "rb"))
uncleaved1=pickle.load(open("uncleaved1.p", "rb"))
uncleaved2=pickle.load(open("uncleaved2.p", "rb"))
uncleaved3=pickle.load(open("uncleaved3.p", "rb"))

HDVtotalnum=0
for key in PercentCleaved1:
    if key in PercentCleaved2 and key in PercentCleaved3:
        HDVtotalnum+=1

HDVtotal={}       
for key in mutnumHDVLIGLib14:
        if key in PercentCleaved1 and key in PercentCleaved2 and key in PercentCleaved3:
            HDVtotal[key]=numpy.mean([PercentCleaved1[key],PercentCleaved2[key],PercentCleaved3[key]])
        else:
            HDVtotal[key]=0
print ("HDV total : ", HDVtotalnum)


normHDVtotal={}
for key in HDVtotal:
    normHDVtotal[key]=(HDVtotal[key]/HDVtotal['CATCGTCCCCGGAC'])*1.25 #normalizes HDV data to HDV7=1.25
pickle.dump(normHDVtotal, open("normHDVtotal.p", "wb"))
#%%
#LIG----------------------------------------------------------------------------------
ribozyme_start = ('GACTCCCA', 'GAATCCCA')
rep=1
while rep!=4:
    ligated={}
    for sequence in open("LIG"+str(rep)+".fastq"):
        sequence=sequence.rstrip()
        if 'AAGCATCTAAGCATCTCAAGCAAACCAGTCG' in sequence and (sequence.startswith('A') or sequence.startswith('C') or sequence.startswith('T') or sequence.startswith('G')):
            if '7' not in sequence and ')' not in sequence and 'N' not in sequence and '<' not in sequence:         
                if 'TCCCATTAG' in sequence and 'GCGGCGGGAGTTG' in sequence and 'AGGGAGGAA' in sequence and 'CCGCCTCCT' in sequence and ('GACTCCCA' in sequence or 'GAATCCCA' in sequence) and 'CTATAGGA' not in sequence:    
                    while not sequence.startswith(ribozyme_start):
                        sequence=sequence.replace(' ', '')[1:]
                    sequence=sequence[2:3]+sequence[12:13]+sequence[17:18]+sequence[27:28]+sequence[41:42]+sequence[44:45]+sequence[54:55]+sequence[57:58]+sequence[61:63]+sequence[67:68]+sequence[72:73]+sequence[74:76]       
                    if (sequence[:1]=='A' or sequence[:1]=='C') and (sequence[1:2]=='A' or sequence[1:2]=='G') and (sequence[2:3]=='G' or sequence[2:3]=='T') and (sequence[3:4]=='G' or sequence[3:4]=='C') and (sequence[4:5]=='G' or sequence[4:5]=='C') and (sequence[5:6]=='G' or sequence[5:6]=='T') and (sequence[6:7]=='G' or sequence[6:7]=='C') and (sequence[7:8]=='C' or sequence[7:8]=='T') and (sequence[8:9]=='C' or sequence[8:9]=='T') and (sequence[9:10]=='C' or sequence[9:10]=='T') and (sequence[10:11]=='G' or sequence[10:11]=='A') and (sequence[11:12]=='G' or sequence[11:12]=='C') and (sequence[12:13]=='A' or sequence[12:13]=='C') and (sequence[13:14]=='G' or sequence[13:14]=='C'):        
                        if sequence in ligated:
                                ligated[sequence]+=1
                        elif sequence not in ligated:
                                ligated[sequence]=1
    print ("Number Ligated "+str(rep), ":", len(ligated))
    pickle.dump(ligated, open("ligated"+str(rep)+".p", "wb"))
    rep+=1

ligated1=pickle.load(open("ligated1.p", "rb"))
ligated2=pickle.load(open("ligated2.p", "rb"))
ligated3=pickle.load(open("ligated3.p", "rb"))

LIGtotalnum=0
for key in ligated1:
    if key in ligated2 and key in ligated3:
        LIGtotalnum+=1

LIGtotal={}       
for key in mutnumHDVLIGLib14:
        if key in ligated1 and key in ligated2 and key in ligated3:
            LIGtotal[key]=numpy.mean([ligated1[key],ligated2[key],ligated3[key]])
        else:
            LIGtotal[key]=0
print ("LIG total : ", LIGtotalnum)
pickle.dump(LIGtotal, open("LIGtotal.p", "wb"))

#%%
#unLIG------------------------------------------------------------------------------------
ribozyme_start=('GACTCCCA', 'GAATCCCA')
rep=1
while rep!=4:
    unligated={}
    for sequence in open("unLIG"+str(rep)+".fastq"):
        sequence=sequence.rstrip()
        if (sequence.startswith('A') or sequence.startswith('C') or sequence.startswith('T') or sequence.startswith('G')):
            if '7' not in sequence and ')' not in sequence and 'N' not in sequence and '<' not in sequence:         
                if 'TCCCATTAG' in sequence and 'GCGGCGGGAGTTG' in sequence and 'AGGGAGGAA' in sequence and 'CCGCCTCCT' in sequence and ('GACTCCCA' in sequence or 'GAATCCCA' in sequence) and 'CTATAGGA' not in sequence:    
                    while not sequence.startswith(ribozyme_start):
                        sequence=sequence.replace(' ', '')[1:]
                    sequence=sequence[2:3]+sequence[12:13]+sequence[17:18]+sequence[27:28]+sequence[41:42]+sequence[44:45]+sequence[54:55]+sequence[57:58]+sequence[61:63]+sequence[67:68]+sequence[72:73]+sequence[74:76]       
                    if (sequence[:1]=='A' or sequence[:1]=='C') and (sequence[1:2]=='A' or sequence[1:2]=='G') and (sequence[2:3]=='G' or sequence[2:3]=='T') and (sequence[3:4]=='G' or sequence[3:4]=='C') and (sequence[4:5]=='G' or sequence[4:5]=='C') and (sequence[5:6]=='G' or sequence[5:6]=='T') and (sequence[6:7]=='G' or sequence[6:7]=='C') and (sequence[7:8]=='C' or sequence[7:8]=='T') and (sequence[8:9]=='C' or sequence[8:9]=='T') and (sequence[9:10]=='C' or sequence[9:10]=='T') and (sequence[10:11]=='G' or sequence[10:11]=='A') and (sequence[11:12]=='G' or sequence[11:12]=='C') and (sequence[12:13]=='A' or sequence[12:13]=='C') and (sequence[13:14]=='G' or sequence[13:14]=='C'):        
                        if sequence in unligated:
                                unligated[sequence]+=1
                        elif sequence not in unligated:
                                unligated[sequence]=1
    print ("Number Unligated "+str(rep), ":", len(unligated))
    pickle.dump(unligated, open("unligated"+str(rep)+".p", "wb"))
    rep+=1

unligated1=pickle.load(open("unligated1.p", "rb"))
unligated2=pickle.load(open("unligated2.p", "rb"))
unligated3=pickle.load(open("unligated3.p", "rb"))


unLIGtotalnum=0
for key in unligated1:
    if key in unligated2 and key in unligated3:
        unLIGtotalnum+=1

unLIGtotal={}       
for key in mutnumHDVLIGLib14:
    unLIGtotal[key]=numpy.mean([unligated1.get(key,0),unligated2.get(key,0),unligated3.get(key,0)])
print ("unLIG total : ", unLIGtotalnum)
pickle.dump(unLIGtotal, open("unLIGtotal.p", "wb"))

LIGmean=numpy.mean(([sum(ligated1.values()),sum(ligated2.values()),sum(ligated3.values())]))
normLIGtotal={}
for key in LIGtotal:
    normLIGtotal[key]=LIGtotal[key]/LIGmean#determines %sequence of total

unLIGmean=numpy.mean([sum(unligated1.values()),sum(unligated2.values()),sum(unligated3.values())])
normunLIGtotal={}    
for key in unLIGtotal:
    normunLIGtotal[key]=unLIGtotal[key]/unLIGmean #determines %sequence of total

LIGunLIG={}
for key in LIGtotal:
    LIGunLIG[key]=(normLIGtotal[key]/normunLIGtotal[key]) #determines %LIG
pickle.dump(LIGunLIG, open("LIGunLIG.p", "wb"))   

normLIGunLIG={}
for key in LIGunLIG:
    normLIGunLIG[key]=LIGunLIG[key]/LIGunLIG['AAGGGGGCTTGCCC'] #normalizes all %LIG to LIG12=1
pickle.dump(normLIGunLIG, open("normLIGunLIG.p", "wb"))

import csv       
with open('HDV-LIGfull.csv','w', newline='') as csvfile:
    for key in mutnumHDVLIGLib14:
        if key in LIGunLIG and key in HDVtotal:
            spamwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow([key, binarysequences[key], cleaved1.get(key,0), uncleaved1.get(key,0), cleaved2.get(key,0), uncleaved2.get(key,0), cleaved3.get(key,0), uncleaved3.get(key,0), HDVtotal[key],normHDVtotal[key], ligated1.get(key,0), ligated2.get(key,0), ligated3.get(key,0), unligated1.get(key,0), unligated2.get(key,0), unligated3.get(key,0), LIGunLIG.get(key,0), normLIGunLIG.get(key,0), mutnumHDVLIGLib14[key]])

with open('HDV-LIGdata.csv','w', newline='') as csvfile:
    for key in mutnumHDVLIGLib14:
        if key in LIGunLIG and key in HDVtotal:
            spamwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow([key, binarysequences[key], normHDVtotal[key], normLIGunLIG.get(key,0), mutnumHDVLIGLib14[key]])