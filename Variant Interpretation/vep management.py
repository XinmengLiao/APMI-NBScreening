#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%Cell 1 Define the path for parameters and databases
import sys
import pandas as pd
import numpy as np
import os
import gzip
import re

print("Newborn vcf files")

run_GWAS = False
run_Pharmaco = False
run_traits = True
run_simplify = False
Upload_To_Cloud = False

fileName = sys.argv[1]
outFile = sys.argv[2]

geneBaseFile = 'Genedb1213.txt'
diseaesDBfile = 'diseaseDB_0605.txt'
OMIM_inheritance_DBfile = 'pheno_OMIM_all.txt'
ontology_file = "genedb.ontology.all0307.csv"

#Read disease and inheritance DB
DiseaseDB = pd.read_csv(diseaesDBfile, sep="\t",encoding="ISO-8859-1")
DiseaseDB= DiseaseDB.replace(np.nan,"No info")
OMIM_Inheritance_DB = pd.read_csv(OMIM_inheritance_DBfile,sep="\t",dtype={"phenotypeMimNumber": str})
OMIM_Inheritance_DB['inheritances'] = OMIM_Inheritance_DB['inheritances'].replace(np.nan,"Inheritance Not provided by OMIM")

#Read omim, hgnc, decipher files
hgnc_file = 'hgnc.txt'
omim_file = 'omim.txt'
decipher_file = 'decipher.txt'
clingen_file = 'Cingen-variants-2024-12-09.txt'

if run_GWAS == True:
    GWAS_dbfile = 'Merged_GWAS_vcf0627.txt'
    GWAS_db =pd.read_csv(GWAS_dbfile, sep="\t")
    GWAS_db= GWAS_db.replace(np.nan,"No info")
    print("GWAS matching will run")

if run_Pharmaco == True:
    Pharma_dbfile = 'AllPGx_annotation1218.txt'
    Pharma_db =pd.read_csv(Pharma_dbfile, sep="\t")
    #Pharma_db= Pharma_db.replace(np.nan,"No info")
    # Read the haplotype variants data
    HaplotypeFile='All_Haplotype_var1218.txt'
    with open(HaplotypeFile, 'r') as file:
        haplotype_var = file.read().strip().splitlines()
    HaplotypeID='All_Haplotype_rsID1218.txt'
    with open(HaplotypeID, 'r') as file:
        haplotype_rsID = file.read().strip().splitlines()
    print("Pharmaco matching will run")

colNames_CSQ =['Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp',
'cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation',
'DISTANCE','STRAND','FLAGS','VARIANT_CLASS','SYMBOL_SOURCE','HGNC_ID','CANONICAL','REFSEQ_MATCH',
'REFSEQ_OFFSET','GIVEN_REF','USED_REF','BAM_EDIT','SOURCE','SIFT','PolyPhen','HGVS_OFFSET','AF',
'gnomADe_AF','gnomADe_AFR_AF','gnomADe_AMR_AF','gnomADe_ASJ_AF','gnomADe_EAS_AF','gnomADe_FIN_AF',
'gnomADe_MID_AF','gnomADe_NFE_AF','gnomADe_REMAINING_AF','gnomADe_SAS_AF','gnomADg_AF','gnomADg_AFR_AF',
'gnomADg_AMI_AF','gnomADg_AMR_AF','gnomADg_ASJ_AF','gnomADg_EAS_AF','gnomADg_FIN_AF','gnomADg_MID_AF',
'gnomADg_NFE_AF','gnomADg_REMAINING_AF','gnomADg_SAS_AF','MAX_AF','MAX_AF_POPS','CLIN_SIG','SOMATIC','PHENO',
'ada_score','rf_score','REVEL','BayesDel_addAF_pred','BayesDel_addAF_rankscore','BayesDel_addAF_score',
'BayesDel_noAF_pred','BayesDel_noAF_rankscore','BayesDel_noAF_score','REVEL_rankscore','REVEL_score','SpliceAI_pred_DP_AG',
'SpliceAI_pred_DP_AL','SpliceAI_pred_DP_DG','SpliceAI_pred_DP_DL','SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG',
'SpliceAI_pred_DS_DL','SpliceAI_pred_SYMBOL','am_class','am_genome','am_pathogenicity','am_protein_variant','am_transcript_id',
'am_uniprot_id','LoF','LoF_filter','LoF_flags','LoF_info','ClinVar','ClinVar_ID','ClinVar_CLNSIG','ClinVar_CLNDN','ClinVar_CLNHGVS',
'ClinVar_CLNSIGINCL','ClinVar_CLNVC','ClinVar_GENEINFO','ClinVar_CLNDISDB','ClinVar_CLNSIGCONF','ClinVar_CLNREVSTAT',
'ClinVar_CLNDNINCL','Database','Database_Type','Database_SZAID']




#%%Cell 2 start the analysis program for the vep annotated vcf files, generated report A B C and D || E for the traits

check_male_flag = False
file = gzip.open(fileName,'rt')
tLine = file.readline()
i = 0
reportA ,reportB,reportC,reportD,reportE = [], [], [], [], []
while tLine.endswith('\n'):
    i = i+1
    iContent = tLine.replace('\n','').split('\t')

    if '#' in iContent[0]:
        if 'ID=CSQ' in tLine:
            annoText = iContent[0].split('Format: ')
            colNames_CSQ = annoText[1].replace('">','')
            colNames_CSQ = colNames_CSQ.split('|')
        elif iContent[0] == '#CHROM':
            headings = iContent
    else:
        iText = [s for s in iContent[headings.index('INFO')].split(';') if 'CSQ=' in s]
        iText = iText[0].replace('CSQ=','').split(',')
        saveFlag1, saveFlag2,saveFlag3,saveFlag4, saveFlag5,saveFlag6 = False, False, False, False, False, False


        for j in range(0,len(iText)):
            jText = iText[j].split('|')

            ichr = iContent[headings.index('#CHROM')].split('_')[0]
            ipos = iContent[headings.index('POS')]
            iref = iContent[headings.index('REF')]
            ialt = iContent[headings.index('ALT')]
            ivariation3 = f"{ipos}_{iref}_{ialt}"
            ivariation4 = f"{ichr}_{ipos}_{iref}_{ialt}"

            is_clinvar = jText[colNames_CSQ.index('Database_SZAID')] != '' and 'ClinP_LP_var' in jText[colNames_CSQ.index('Database_Type')].split('&')

            if is_clinvar:
                saveFlag1 = 'ClinP_LP_var'
                
            if jText[colNames_CSQ.index('MAX_AF')] == '' : 
                jText[colNames_CSQ.index('MAX_AF')] = 0
            jText[colNames_CSQ.index('MAX_AF')] = float(jText[colNames_CSQ.index('MAX_AF')])
            if not saveFlag1 and not saveFlag6 and jText[colNames_CSQ.index('MAX_AF')] < 0.05:
                ##Update: to clear the logics here need to add brackets
                if (jText[colNames_CSQ.index('IMPACT')] == 'HIGH' or (
                        jText[colNames_CSQ.index('ada_score')] != '' and float(jText[colNames_CSQ.index('ada_score')])>0.6) or (
                        jText[colNames_CSQ.index('rf_score')] != '' and float(jText[colNames_CSQ.index('rf_score')])>0.6) or (
                        jText[colNames_CSQ.index('REVEL')]!= '' and float(jText[colNames_CSQ.index('REVEL')])>0.75) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_AL')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_AL')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_DG')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_DG')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_DL')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_DL')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_AG')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_AG')])>0.5) or (
                        jText[colNames_CSQ.index('BayesDel_addAF_score')]!= '' and float(jText[colNames_CSQ.index('BayesDel_addAF_score')])>0.0692655) or (
                        jText[colNames_CSQ.index('BayesDel_noAF_score')]!= '' and float(jText[colNames_CSQ.index('BayesDel_noAF_score')])>-0.0570105) or (
                        jText[colNames_CSQ.index('am_class')]== 'likely_pathogenic' and float(jText[colNames_CSQ.index('am_pathogenicity')])> 0.564) or (
                        jText[colNames_CSQ.index('LoF')]== 'HC')):
                    saveFlag2 ="Putative_var"
            
            if jText[colNames_CSQ.index('Database_SZAID')] != '' and 'GWAS_var' in jText[colNames_CSQ.index('Database_Type')].split('&'):
                saveFlag3 = "GWAS_var"
            if jText[colNames_CSQ.index('Database_SZAID')] != '' and 'Pharma_var' in jText[colNames_CSQ.index('Database_Type')].split('&'):
                saveFlag4 = "Pharma_var"

            if run_Pharmaco == True:  
                if ivariation4 in haplotype_var:
                    saveFlag4 = "Pharma_var"    

                rsID = jText[colNames_CSQ.index('Existing_variation')]
                rsID = rsID.split('&')
                if True in [part in haplotype_rsID for part in rsID]:
                    saveFlag4 = "Pharma_var" 
            
            # Check gender 
            if 'chrY' in iContent[headings.index('#CHROM')] :
                check_male_flag = True

        if saveFlag1:
            reportA.append(tLine)
        if saveFlag2:
            reportB.append(tLine)
        if saveFlag3:
            reportC.append(tLine)
        if saveFlag4:
            reportD.append(tLine)
        if saveFlag5:
            reportE.append(tLine)

    if i%100000 == 0:
        print(str(i)+' lines processed!')
    tLine = file.readline()
    
variant_count = i
gender = "Male" if check_male_flag else "Female"

if gender == "Male" :
    print("Male")
else:
    print("Female")

file.close()
print('Cell 2 VEP annotated File processing done! Now start to map GeneDB and DiseaseDB')



#%%Cell 3 read gene db and disease db

print("running cell 3,reading gene db and disease db")
# read gene based database
ConfidenceLevel = []
TargetGroup = []
geneBasedRef = []
Disease = []
SZAdiseaseID_GDB = []
inh = []
mark = []
mim = []

with open(geneBaseFile,'r') as f:
    for line in f:
        line = line.replace('\n','')
        temp = line.split('\t')
        ConfidenceLevel.append(temp[0])
        TargetGroup.append(temp[1])
        geneBasedRef.append(temp[2])
        Disease.append(temp[3])
        SZAdiseaseID_GDB.append(temp[4])
        inh.append(temp[5])
        mark.append(temp[6])
        mim.append(temp[7])
geneBasedRefHeadings = [ConfidenceLevel[0],TargetGroup[0],geneBasedRef[0],Disease[0],inh[0],mark[0],mim[0]]

# Remove the first line for each col, which is the heading 
geneBasedRef.pop(0)
Disease.pop(0)
ConfidenceLevel.pop(0)
TargetGroup.pop(0)
SZAdiseaseID_GDB.pop(0)
inh.pop(0)
mim.pop(0)
mark.pop(0)

# read disease based database
DiseaseID_DSDB = []
DiseaseName_DSDB = []
count=0
with open(diseaesDBfile,'r',encoding = "ISO-8859-1") as f:
    for line in f:
        count+=1
        line = line.replace('\n','')
        temp = line.split('\t')
        DiseaseID_DSDB.append(temp[0])
        DiseaseName_DSDB.append(temp[1])
        #print(line)
DiseaseID_DSDB.pop(0)
DiseaseName_DSDB.pop(0)

newContentHeadings = ['SZAID','Database','Database_type','SZAdiseaseID','ClinVar_CLNSIG','ClinVar_CLNDN','SZAreportCategory','ClinVar_ReviewStar','ClinVar_ID','selCSQ','Zygosity','Genotype']+geneBasedRefHeadings+headings;
newContent = '\t'.join(newContentHeadings)+'\n'


#%%Cell 4 report A is for Clinvar Variants
# Part A

for i in range(0,len(reportA)):
    iText = [s for s in reportA[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
    iREF = reportA[i].split('\t')[3].split(',')
    iALT = reportA[i].split('\t')[4].split(',')
    iGenoTypeList = iREF+iALT
    if "." not in reportA[i].split('\t')[9].split(':')[0]:
        iGenotypInd1 = int(reportA[i].split('\t')[9][0])
        iGenotypInd2 = int(reportA[i].split('\t')[9][2])
        if iGenotypInd2 <= len(iGenoTypeList)-1 :
            iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
        else:
            iGenotype = iREF[0]+'>'+iGenoTypeList[iGenotypInd1]
        iTemp = reportA[i].replace('\n','').split('\t')
        iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
        if iGenotypInd1 == iGenotypInd2:
            iZygo = 'Homozygous'
        elif iGenotypInd1 == 0 or iGenotypInd2 == 0 or iGenotypInd2 > len(iGenoTypeList)-1:
            iZygo = 'Heterozygous'
        if check_male_flag == "Male":
            if "chrX" in reportA[i] and iZygo == 'Heterozygous':
                iZygo = "Hemizygous"
        for j in range(0,len(iText)):
            jText = iText[j].split('|')
            jGenes1 = jText[3]
            jGenes2Temp = jText[colNames_CSQ.index("ClinVar_GENEINFO")].split('&')
            jGenes2 = []
            for k in range(0,len(jGenes2Temp)):
                kTemp = jGenes2Temp[k].split(':')
                jGenes2.append(kTemp[0])
            jGenes = list(set([jGenes1]+jGenes2))
            jGenes = [gene for gene in jGenes if gene != ""]
            
            #For CLN Review Star 
            Review = jText[colNames_CSQ.index('ClinVar_CLNREVSTAT')]
            if Review == 'criteria_provided&_conflicting_classifications':
                ReviewStar = 1
            elif Review == 'criteria_provided&_single_submitter':
                ReviewStar = 1
            elif Review == 'criteria_provided&_multiple_submitters&_no_conflicts':
                ReviewStar = 2
            elif Review == 'reviewed_by_expert_panel':
                ReviewStar = 3
            elif Review == 'practice_guideline':
                ReviewStar = 4
            else:
                ReviewStar = 0
            
            jCLINVARind = jText[colNames_CSQ.index('Database_Type')].split('&').index('ClinP_LP_var') 
            jDatabaseType = jText[colNames_CSQ.index('Database_Type')].split('&')[jCLINVARind]
            jCLNID = jText[colNames_CSQ.index('Database')].split('&')[jCLINVARind]
            jSZAvarID = jText[colNames_CSQ.index('Database_SZAID')].split('&')[jCLINVARind]
            jCLNSIG = jText[colNames_CSQ.index('ClinVar_CLNSIG')]
            jClinVar_CLNDN =  jText[colNames_CSQ.index('ClinVar_CLNDN')].replace('&_',',_').replace('&','|').replace('_',' ') 
            if 'Pathogenic' in jCLNSIG:
                scoreFlag = 5
            elif 'Likely_pathogenic' in jCLNSIG:
                scoreFlag = 4
            elif 'Conflicting_interpretations_of_pathogenicity' in jCLNSIG:
                scoreFlag = 1
            elif 'Conflicting_classifications_of_pathogenicity' in jCLNSIG:
                pathogenicity = jText[colNames_CSQ.index('ClinVar_CLNSIGCONF')].split('&')
                pathogenicity = [re.sub(r'\(\d+\)', '', s) for s in pathogenicity]
                if 'Pathogenic' in pathogenicity or 'Likely_pathogenic' in pathogenicity:    
                    scoreFlag = 2
                else:
                    scoreFlag = 1
            else:
                scoreFlag = 0

            ## 1. gene is in GeneDB
            if scoreFlag > 0 and any([x for x in jGenes if x in geneBasedRef]): 
                jGeneList = [x for x in jGenes if x in geneBasedRef] # store all associated genes as a list 
                jInds = []

                for gene in range(0,len(jGeneList)):
                    kTemp = jGeneList[gene]
                    jInds = jInds +[ii for ii, x in enumerate(geneBasedRef) if x == kTemp] # Find gene location in GeneDB


                # Extract the ClinVar disease in each line 
                jClinVar_CLNDN =  jText[colNames_CSQ.index('ClinVar_CLNDN')].replace('&_',',_').replace('&','|').replace('_',' ')
                kClinVarDiseases = jClinVar_CLNDN.split('|') # All the ClinVar disease in each line 
                jVarSZAdiseaseIDs = []
                jVarDiseaseNames = []
                for c in range(0,len(kClinVarDiseases)):
                    kClinVarDisease_temp = kClinVarDiseases[c].replace(',','')
                    if len([ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]) == 1: 
                        kIndex = [ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp == x.replace(',','')]
                        for element in kIndex:
                            jVarSZAdiseaseIDs.append(DiseaseID_DSDB[element])
                            jVarDiseaseNames.append(DiseaseName_DSDB[element])

                for k in range(0,len(jInds)):
                    if SZAdiseaseID_GDB[jInds[k]] in jVarSZAdiseaseIDs: # if the GeneDB disease is in DiseaseDB  disease (exactly match)
                        if ConfidenceLevel[jInds[k]] == 'High_confidence':
                            SZAscore = 10 + scoreFlag
                        elif ConfidenceLevel[jInds[k]] == 'Moderate_confidence':
                            SZAscore = 5 + scoreFlag
                        elif ConfidenceLevel[jInds[k]] == 'Low_confidence':
                            SZAscore = scoreFlag
                        else:
                            sys.exit('Unknown confidence score level for the report!')
                        jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]
                        
                        jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                        jCLNID = jText[colNames_CSQ.index('Database')]
                        newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                                    Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]] +iTemp)+'\n'
                                        
                    else:
                        GeneDB_szaID = SZAdiseaseID_GDB[jInds[k]]
                        GeneDB_Disease = Disease[jInds[k]]
                        GeneDB_Disease_word = set(GeneDB_Disease.upper().split())

                        if len(jVarDiseaseNames) > 0:
                            for name in range(len(jVarDiseaseNames)):
                                clinvar_disease_word = set(jVarDiseaseNames[name].upper().split())
                                matches = False

                                if clinvar_disease_word.issubset(GeneDB_Disease_word):
                                    matches = True
                                elif GeneDB_Disease_word.issubset(clinvar_disease_word):
                                    matches = True
                            
                                if matches == True:
                                    if ConfidenceLevel[jInds[k]] == 'High_confidence':
                                    #print('SZAdiseaseID_GDB[jInds[k]]=',SZAdiseaseID_GDB[jInds[k]])
                                    #print(' jVarSZAdiseaseIDs=',jVarSZAdiseaseIDs)
                                        SZAscore = 10 + scoreFlag
                                    elif ConfidenceLevel[jInds[k]] == 'Moderate_confidence':
                                        SZAscore = 5 + scoreFlag
                                    elif ConfidenceLevel[jInds[k]] == 'Low_confidence':
                                        SZAscore = scoreFlag
                                    else:
                                        sys.exit('Unknown confidence score level for the report!')
                                    jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]

                                    jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                                    jCLNID = jText[colNames_CSQ.index('Database')]                                    
                                    newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                                                Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]] +iTemp)+'\n'

                                else:
                                ########here, if SZAdiseaseID_GDB[jInds[k]] in jVarSZAdiseaseIDs is false, then no SZAscore value
                                    SZAscore = scoreFlag
                                    jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]

                                    jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                                    jCLNID = jText[colNames_CSQ.index('Database')]
                                    newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar),jCLNID,iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                                                    Disease[jInds[k]],inh[jInds[k]],mark[jInds[k]], mim[jInds[k]]]+ iTemp)+'\n'
                                   
                        else: # for variants do not have ClinVar disease
                            SZAscore = scoreFlag 
                            newContent = newContent+'\t'.join(['NA',jCLNID,jDatabaseType,SZAdiseaseID_GDB[jInds[k]],jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID, iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                            Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]] +iTemp)+'\n' 


                            
            elif scoreFlag >0 and any(jGenes): # If there is jGenes, but not included in GeneDB
                #sys.exit('Found a gene set with '+', '.join(jGenes)+' that is not included in GeneDB ')
                SZAscore = scoreFlag
                jClinVar_CLNDN =  jText[colNames_CSQ.index('ClinVar_CLNDN')].replace('&_',',_').replace('&','|').replace('_',' ')
                kClinVarDiseases = jClinVar_CLNDN.split('|')
                #print(kClinVarDiseases)
                jVarSZAdiseaseIDs = []
                jVarDiseaseNames = []
                jSZADiseaseID = []
                jDisease = []
                for k in range(0,len(kClinVarDiseases)):
                    #kClinVarDisease_temp = kClinVarDiseases[k]
                    kClinVarDisease_temp = kClinVarDiseases[k].replace(",","")
                    if len([ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]) == 1:
                        kIndex = [ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]
                        jSZADiseaseID_tmp = DiseaseID_DSDB[kIndex[0]]
                        jSZADiseaseID.append(jSZADiseaseID_tmp)
                        jDisease_tmp = DiseaseName_DSDB[kIndex[0]]
                        jDisease.append(jDisease_tmp)
                    else:
                        jSZADiseaseID_tmp = ''
                        jSZADiseaseID.extend(jSZADiseaseID_tmp)
                        jDisease_tmp = ''
                        jDisease.append(jDisease_tmp)
                jDisease = [x for x in jDisease if x != '']
                jSZADiseaseID = [x for x in jSZADiseaseID if x != '']
                if len(jDisease) == 1: 
                    jDisease = ''.join(jDisease)
                    jSZADiseaseID = ''.join(jSZADiseaseID)

                    jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                    jCLNID = jText[colNames_CSQ.index('Database')]
                    newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),jDisease,"None","None","None"]+iTemp)+'\n'
                   
                
                elif len(jDisease) > 1:
                    for m in range(len(jDisease)):
                        mjDisease = ''.join(jDisease[m])
                        mjSZADiseaseID = ''.join(jSZADiseaseID[m])
                        newContent = newContent+'\t'.join(['NA','NA','HGMD',mjSZADiseaseID,jCLNSIG,'NA',str(0),str(ReviewStar), 'NA',iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),mjDisease,"None","None","None"]+iTemp)+'\n'
                else: # len(jDisease) == 0, for HGMD variant without ClinVar disease 
                    newContent = newContent+'\t'.join(['NA',jCLNID,jDatabaseType,"NA",jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),"None","None","None","None"]+iTemp)+'\n'

            elif scoreFlag >0:
                SZAscore = scoreFlag
                jClinVar_CLNDN =  jText[colNames_CSQ.index('ClinVar_CLNDN')].replace('&_',',_').replace('&','|').replace('_',' ')
                kClinVarDiseases = jClinVar_CLNDN.split('|')
                jVarSZAdiseaseIDs = []
                jVarDiseaseNames = []
                jSZADiseaseID = []
                jDisease = []
                for k in range(0,len(kClinVarDiseases)):
                    kClinVarDisease_temp = kClinVarDiseases[k].replace(",","")
                    if len([ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]) == 1:
                        kIndex = [ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]

                        jSZADiseaseID_tmp = DiseaseID_DSDB[kIndex[0]]
                        jSZADiseaseID.append(jSZADiseaseID_tmp)
                        jDisease_tmp = DiseaseName_DSDB[kIndex[0]]
                        jDisease.append(jDisease_tmp)
                    else:
                        jSZADiseaseID_tmp = ''
                        jSZADiseaseID.extend(jSZADiseaseID_tmp)
                        jDisease_tmp = ''
                        jDisease.append(jDisease_tmp)
                jDisease = [x for x in jDisease if x != '']
                jSZADiseaseID = [x for x in jSZADiseaseID if x != '']
                if len(jDisease) == 1: 
                    jDisease = ''.join(jDisease)
                    jSZADiseaseID = ''.join(jSZADiseaseID)
                    jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                    jCLNID = jText[colNames_CSQ.index('Database')]
                    newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),jDisease,"None","None","None"]+iTemp)+'\n'
                
                elif len(jDisease) > 1:
                    for m in range(len(jDisease)):
                        mjDisease = ''.join(jDisease[m])
                        mjSZADiseaseID = ''.join(jSZADiseaseID[m])
                        jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                        jCLNID = jText[colNames_CSQ.index('Database')]
                        newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,mjSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),mjDisease,"None","None","None"]+iTemp)+'\n'
                else:
                    newContent = newContent+'\t'.join(['NA',jCLNID,jDatabaseType,"NA",jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),"None","None","None","None"]+iTemp)+'\n'


## Compound Heterozygous analyses (ClinVar and HGMD variants will both count)
clinvar_headings = newContent.split('\n')[0].split('\t')
line = [clinvar.split('\t') for clinvar in newContent.split('\n')[1:]]
clinvar_df = pd.DataFrame(data=line, columns=clinvar_headings)
clinvar_df["Compound_heterozygous"] = "No"
compound_tmp_df = clinvar_df[['Genes','#CHROM','POS','REF','ALT',clinvar_df.columns[28]]].drop_duplicates().dropna()
compound_tmp_df.iloc[:, 5] = compound_tmp_df.iloc[:, 5].apply(lambda x: x.split(':')[0])
compound_tmp_df['Genotype'] = None
compound_tmp_df = compound_tmp_df.reset_index(drop=True)

#Extract the variation and zygosity for compound heterozygous 
for i in range(len(compound_tmp_df)):
    genotype0 = compound_tmp_df.iloc[i,5][0]
    if len(compound_tmp_df.iloc[i,5]) > 1:
        genotype1 = compound_tmp_df.iloc[i,5][1]
        genotype2 = compound_tmp_df.iloc[i,5][2]
    else:
        genotype1 = genotype0
        genotype2 = genotype0
    if genotype0 != genotype2:
        if genotype1 == "/":
            compound_tmp_df.at[i,"Genotype"] = "Heterozygous(Unphased)"
        elif genotype1 == "|":
            compound_tmp_df.at[i,"Genotype"]  = "Heterozygous(Phased)"

compound_tmp_df = compound_tmp_df.dropna()
grouped = compound_tmp_df.groupby('Genes').filter(lambda x: len(x) > 1)

unique_genes = grouped['Genes'].unique()
for gene in unique_genes:
    tmp_df = grouped[grouped['Genes'] == gene]
    genotype = tmp_df[tmp_df['Genes'] == gene]['Genotype']
    if 'Heterozygous(Unphased)' in genotype.values:
        for index, row in tmp_df.iterrows():
            ref_tmp = row["REF"]
            chrom_tmp = row['#CHROM']
            alt_tmp = row['ALT']
            pos_tmp = row['POS']
            index = clinvar_df[(clinvar_df['Genes'] == gene) & 
                            (clinvar_df['REF'] == ref_tmp) &
                            (clinvar_df['#CHROM'] == chrom_tmp) &
                            (clinvar_df['POS'] == pos_tmp) &
                            (clinvar_df['ALT'] == alt_tmp)].index
            if len(index) > 0:
                clinvar_df.loc[index, 'Compound_heterozygous'] = 'Uncertained ' + row.iloc[5]
    else:
        clinvar_df.loc[index, 'Compound_heterozygous'] = 'Yes'


## Change df into NewContent
newContent = '\t'.join(newContentHeadings + ["Compound_heterozygous"])+'\n'
for index, row in clinvar_df.iterrows():
    newContent += '\t'.join(row.astype(str)) + '\n'
                        

#%%Cell 5 report C GWAS

if run_GWAS== True:
    from datetime import datetime
    newLineHeadings = ['Trait','Genotype','Zygosity','risk_genotype','SZAvarID','GWASID','OR','Pval','PUBMED','Study accession ID','MAPPED_GENE',"CONTEXT"];
    newLine = '\t'.join(newLineHeadings)+'\n'
    with open(outFile.replace('.txt','_GWAS.txt'),'w') as f:
        f.write(newLine)
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    print('processing ' + str(len(reportC)) + ' lines of GWAS variants' )
    # GWAS part
    for i in range(0,len(reportC)):
        current_time = now.strftime("%H:%M:%S")
        iText = [s for s in reportC[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
        iREF = reportC[i].split('\t')[3].split(',')
        iALT = reportC[i].split('\t')[4].split(',')
        iGenoTypeList = iREF+iALT

        if "." not in reportC[i].split('\t')[9].split(':')[0]:
            iGenotypInd1 = int(reportC[i].split('\t')[9][0])
            iGenotypInd2 = int(reportC[i].split('\t')[9][2])
            try:
                iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
            except:
                iGenotype = iREF[0]+'>'+iALT[0]+':'+reportC[i].split('\t')[9][0]+reportC[i].split('\t')[9][1]+reportC[i].split('\t')[9][2]
            iTemp = reportC[i].replace('\n','').split('\t')
            iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
            if iGenotypInd1 == iGenotypInd2:
                iZygo = 'Homozygous'
            elif iGenotypInd1 == 0 or iGenotypInd2 == 0:
                iZygo = 'Heterozygous'
            else:
                iZygo = 'Compound heterozygous'
            jText = iText[0].split('|')
            if jText[colNames_CSQ.index('Database_Type')]!= '':
                jGWASind = jText[colNames_CSQ.index('Database_Type')].split('&').index('GWAS_var')
                jDatabaseType = jText[colNames_CSQ.index('Database_Type')].split('&')[jGWASind]
                jGWASID = jText[colNames_CSQ.index('Database')].split('&')[jGWASind]
                jSZAvarID = jText[colNames_CSQ.index('Database_SZAID')].split('&')[jGWASind]
                GWAS_db_filtered = GWAS_db[GWAS_db['SZAID'] == jSZAvarID]
                for k, row in GWAS_db_filtered.iterrows():
                    jOR = row['OR.or.BETA']
                    jPval = row['P.VALUE']
                    jrisk_genotype = row['risk_genotype']
                    jPUBMED = row['PUBMEDID']
                    jTrait = row['DISEASE.TRAIT']
                    jMappedGene = row['MAPPED_GENE']
                    jVarType = row['CONTEXT']
                    jStudy_accession = row['STUDY.ACCESSION']
                    if iGenotype == jrisk_genotype and jPval <= 5e-8:
                        newLine_temp = '\t'.join([jTrait, iGenotype, iZygo, jrisk_genotype, jSZAvarID, jGWASID, str(jOR), str(jPval), str(jPUBMED), str(jStudy_accession), jMappedGene, jVarType]) + '\n'
                        with open(outFile.replace('.txt', '_GWAS.txt'), 'a') as f:
                            f.write(newLine_temp)

        if i%1000 == 0:
            print(str(i)+' lines processed!')

    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)


#%% Cell 6 (new) report D for pharmaco-annotation

if run_Pharmaco == True:

    from datetime import datetime
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    print('processing ' + str(len(reportD)) + ' lines of pharmaco variants' )

    for i in range(0,len(reportD)):
        current_time = now.strftime("%H:%M:%S")
        iText = [s for s in reportD[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
        iREF = reportD[i].split('\t')[3].split(',')
        iALT = reportD[i].split('\t')[4].split(',')
        iGenoTypeList = iREF+iALT

        if "." not in reportD[i].split('\t')[9].split(':')[0]:
            iGenotypInd1 = int(reportD[i].split('\t')[9][0])
            iGenotypInd2 = int(reportD[i].split('\t')[9][2])
            genotype1 = iGenoTypeList[iGenotypInd1] + iGenoTypeList[iGenotypInd2]
            genotype_sorted = ''.join(sorted(genotype1))
            try:
                iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
            except:
                iGenotype = iREF[0]+'>'+iALT[0]+':'+reportD[i].split('\t')[9][0]+reportD[i].split('\t')[9][1]+reportD[i].split('\t')[9][2]
            iTemp = reportD[i].replace('\n','').split('\t')
            iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
            if iGenotypInd1 == iGenotypInd2:
                iZygo = 'Homozygous'
            elif iGenotypInd1 == 0 or iGenotypInd2 == 0:
                iZygo = 'Heterozygous'
            else:
                iZygo = 'Compound heterozygous'
            jText = iText[0].split('|')
            if jText[colNames_CSQ.index('Database_Type')]!= '':
                content = jText[colNames_CSQ.index('Database_Type')].split('&')
                if 'Pharma_var' in content:
                    jPharmaind = jText[colNames_CSQ.index('Database_Type')].split('&').index('Pharma_var')
                    jDatabaseType = jText[colNames_CSQ.index('Database_Type')].split('&')[jPharmaind]
                    jPharmaID = jText[colNames_CSQ.index('Database')].split('&')[jPharmaind]
                else:
                    IDs = jText[colNames_CSQ.index('Existing_variation')].split('&')
                    index = next((i for i, ID in enumerate(IDs) if ID.startswith('rs')), None)
                    if index != None:
                        jPharmaID = IDs[index]
                    else:
                        jPharmaID = 'na'
            else:
                IDs = jText[colNames_CSQ.index('Existing_variation')].split('&')
                if IDs != ['']:
                    index = next((i for i, ID in enumerate(IDs) if ID.startswith('rs')), None)
                    if index != None:
                        jPharmaID = IDs[index]
                    else: 
                        jPharmaID = 'na'
                else:
                    jPharmaID = 'na'
                
            jPharmaID = str(jPharmaID)


            ichr = iTemp[0]
            ipos = iTemp[1]
            iref = iTemp[3]
            ialt = iTemp[4]
            ivariation = f"{ichr}_{ipos}_{iref}_{ialt}"


            for r in range(len(Pharma_db)):
                tmp_id = Pharma_db.loc[Pharma_db.index[r], 'rsID']
                if jPharmaID == tmp_id:
                    if '*' in Pharma_db.loc[Pharma_db.index[r], 'Genotype.Allele']:
                        Pharma_db.loc[Pharma_db.index[r],"Genotype"] = iGenoTypeList[iGenotypInd1] + '/' + iGenoTypeList[iGenotypInd2]
                        Pharma_db.loc[Pharma_db.index[r],"Zygous"]= iZygo
                    else:
                        Pharma_db.loc[Pharma_db.index[r],"Genotype"] = genotype_sorted
                        Pharma_db.loc[Pharma_db.index[r],"Zygous"]= iZygo
                if ivariation == tmp_id:
                    if '*' in Pharma_db.loc[Pharma_db.index[r], 'Genotype.Allele']:
                        Pharma_db.loc[Pharma_db.index[r],"Genotype"] = iGenoTypeList[iGenotypInd1] + '/' + iGenoTypeList[iGenotypInd2]
                        Pharma_db.loc[Pharma_db.index[r],"Zygous"]= iZygo
                    else:
                        Pharma_db.loc[Pharma_db.index[r],"Genotype"] = genotype_sorted
                        Pharma_db.loc[Pharma_db.index[r],"Zygous"]= iZygo

    def filter_genotypes(group):
        return not (group['Genotype'] == 'a').any()
    Pharma_db_final = Pharma_db.groupby('Variant.Haplotypes').filter(filter_genotypes)

    Pharma_db_final = Pharma_db_final[
    ~((Pharma_db_final['Variant.Haplotypes'].str.startswith('rs')) & 
      (Pharma_db_final['Genotype'] != Pharma_db_final['Genotype.Allele']))]
                
    pgx_filename =outFile.replace('.txt','_PGx.txt') 

    Pharma_db_final.to_csv(pgx_filename, index=False,sep='\t')
  
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)

#%% Cell 7 Part B processing and generate outFile
varCount = 0
for i in range(0,len(reportB)):
    iText = [s for s in reportB[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
    iREF = reportB[i].split('\t')[3].split(',')
    iALT = reportB[i].split('\t')[4].split(',')
    iGenoTypeList = iREF+iALT

    if "." not in reportB[i].split('\t')[9].split(':')[0]:
        iGenotypInd1 = int(reportB[i].split('\t')[9][0])
        iGenotypInd2 = int(reportB[i].split('\t')[9][2])
        try:
            iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
        except:
            iGenotype = iREF[0]+'>'+iALT[0]+':'+reportB[i].split('\t')[9][0]+reportB[i].split('\t')[9][1]+reportB[i].split('\t')[9][2]
        if iGenotypInd1 == iGenotypInd2:
            iZygo = 'Homozygous'
        elif iGenotypInd1 == 0 or iGenotypInd2 == 0:
            iZygo = 'Heterozygous'
        else:
            iZygo = 'Compound heterozygous'
        if check_male_flag == "Male":
            if "chrX" in reportB[i] and iZygo == 'Heterozygous':
                iZygo = "Hemizygous"
        varFlag = 0
        iTemp = reportB[i].replace('\n','').split('\t')
        iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
        for j in range(0,len(iText)):
            jText = iText[j].split('|')
            if (jText[colNames_CSQ.index('MAX_AF')] == '' or float(jText[colNames_CSQ.index('MAX_AF')])<0.05):
                   
                if (jText[colNames_CSQ.index('IMPACT')] == 'HIGH' or (
                        jText[colNames_CSQ.index('ada_score')] != '' and float(jText[colNames_CSQ.index('ada_score')])>0.6) or (
                        jText[colNames_CSQ.index('rf_score')] != '' and float(jText[colNames_CSQ.index('rf_score')])>0.6) or (
                        jText[colNames_CSQ.index('REVEL')]!= '' and float(jText[colNames_CSQ.index('REVEL')])>0.75) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_AL')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_AL')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_DG')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_DG')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_DL')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_DL')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_AG')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_AG')])>0.5) or (
                        jText[colNames_CSQ.index('BayesDel_addAF_score')]!= '' and float(jText[colNames_CSQ.index('BayesDel_addAF_score')])>0.0692655) or (
                        jText[colNames_CSQ.index('BayesDel_noAF_score')]!= '' and float(jText[colNames_CSQ.index('BayesDel_noAF_score')])>-0.0570105) or (
                        jText[colNames_CSQ.index('am_class')]== 'likely_pathogenic' and float(jText[colNames_CSQ.index('am_pathogenicity')])> 0.564) or (
                        jText[colNames_CSQ.index('LoF')]== 'HC')):

                    jGenes = jText[3].split('&')
                    if any([x for x in jGenes if x in geneBasedRef]):
                        varFlag = 1
                        jGeneList = [x for x in jGenes if x in geneBasedRef]
                        jInds = []
                        for k in range(0,len(jGeneList)):
                            kTemp = jGeneList[k]
                            jInds = jInds +[i for i, x in enumerate(geneBasedRef) if x == kTemp]
                        for k in range(0,len(jInds)):
                            scoreFlag = 2+(jText[colNames_CSQ.index('IMPACT')]=='HIGH')*1
                            if ConfidenceLevel[jInds[k]] == 'High_confidence':
                                SZAscore = 10+scoreFlag
                            elif ConfidenceLevel[jInds[k]] == 'Moderate_confidence':
                                SZAscore = 5+scoreFlag
                            elif ConfidenceLevel[jInds[k]] == 'Low_confidence':
                                SZAscore = scoreFlag
                            else:
                                print(jGenes)
                                sys.exit('Unknown confidence score level for the report!')
                            jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]
                            newContent = newContent+'\t'.join(
                                ['NovelTrans'+headings[9]+'_'+str(varCount+1),'NA','NA',
                                jSZADiseaseID,'NA','NA',str(SZAscore),'NA','NA',iText[j],iZygo,
                                iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],
                                geneBasedRef[jInds[k]],Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]]+iTemp + ["No"])+'\n'
                
                
    if varFlag == 1:
        varCount = varCount+1

print("putative variants matched the GeneDB:",varCount)

# write final report
outfid = open(outFile,'w')
outfid.write(newContent)
outfid.close()

print("outFile generated. ")


outFile_sp = outFile.replace('.txt','_sp.txt')
final_report_sp  = pd.read_csv(outFile, sep="\t")
split_df = final_report_sp['selCSQ'].str.split("|", expand=True)
split_df.columns = colNames_CSQ
final_report_sp = pd.concat([final_report_sp, split_df], axis=1)
final_report_sp.to_csv(outFile_sp, sep = "\t",index = None, header=True,)
outFile_sp_Inheritance = outFile.replace('.txt','_sp_Inheritance.txt')
NewHeadings = []
NewHeadings = list(final_report_sp.columns) + ['Inheritances'] + ['DiseaseInfo']

with open((outFile_sp_Inheritance),'w') as f2:
    f2.write('\t'.join(NewHeadings)+'\n')

#matching the inheritance
with open(outFile_sp,'r') as f1:
    next(f1)
    for line in f1:
        line = line.replace('\n','')
        temp = line.split('\t')
        iSZAdiseaseID = temp[3]
        if 'SZA' not in temp[3]:
            continue
        iGene = temp[31] 
        iOMIM_num = DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['DiseaseMIM'].to_list()[0]
        if iOMIM_num:
            iheritance = OMIM_Inheritance_DB[(OMIM_Inheritance_DB['phenotypeMimNumber']== iOMIM_num) &
                                        (OMIM_Inheritance_DB['approvedGeneSymbol']== iGene)]['inheritances'].to_list()
            if iheritance:
                iheritance = iheritance[0]
            else:
                iheritance = "No matched diseaseMIM in OMIM/OMIM not provided Inheritance"
        iDiseaseInfo = DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['OMIM_Description'].to_list()[0]+ '|' + "Reference:"+ "DiseaseName:" + DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['DiseaseName'].to_list()[0]+ '|' + "DiseaseSource:" + DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['SourceName'].to_list()[0]+'|' +"Disease SourceID:" +DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['SourceID'].to_list()[0]+'|' + "OMIM number:"+ DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['DiseaseMIM'].to_list()[0]
        line = line + "\t" + iheritance +"\t"+iDiseaseInfo
        with open((outFile_sp_Inheritance),'a') as f2:
            f2.write(line+ '\n')

final_report_sp_Inheritance  = pd.read_csv(outFile_sp_Inheritance, sep="\t",dtype=str)
final_report_general_4 = final_report_sp_Inheritance
final_report_general_4.to_csv(outFile.replace('.txt','_sp_Inheritance_4.txt'), sep = "\t",index = None, header=True)

print("the file named _sp_Inheritance generated: XXX_sp_Inheritance_1.txt,XXX_sp_Inheritance_2.txt,XXX_sp_Inheritance_3.txt,XXX_sp_Inheritance_4")


#%% Cell 8 def append_data function
colNames = list(final_report_sp_Inheritance.columns)
def append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID,iinh, imark, imim, iExisting_variation, iCHROM, iPOS,
                                          iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                          iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory,istar,icompound_het,iSymbol):
    line_list.append(iline)
    TargetGroup.append(iTargetGroup)
    Disease.append(iDisease)
    Gene.append(iGene)
    SZAID.append(iSZAID)
    inh.append(iinh)
    mark.append(imark)
    mim.append(imim)
    Existing_variation.append(iExisting_variation)
    CHROM.append(iCHROM)
    POS.append(iPOS)
    Zygosity.append(iZygosity)
    REF.append(iREF)
    ALT.append(iALT)
    ID.append(iID)
    ClinVar_CLNSIG.append(iClinvar_CLNSIG)
    VARIANT_CLASS.append(iVARIANT_CLASS)
    Consequence.append(iConsequence)
    ClinVar_CLNHGVS.append(iClinVar_CLNHGVS)
    Inheritances.append(iheritance)
    VariantINFO.append(iVariantINFO)
    DiseaseInfo.append(iDiseaseInfo)
    SZAreportCategory.append(iSZAreportCategory)
    ReviewStar.append(istar)
    Compound_heterozygous.append(icompound_het)
    Symbol.append(iSymbol)


# create lists
TargetGroup = []
Disease = []
Gene = []
Symbol = []
SZAID = []
inh = []
mark = []
mim = []
Existing_variation =[]
VARIANT_CLASS =[]
Consequence =[]
ClinVar_CLNHGVS=[]
#Genotype = []
CHROM = []
POS = []
ID = []
REF=[]
ALT=[]
Zygosity = []
ConfidenceLevel = []
#SZAdiseaseID = []
ClinVar_CLNSIG=[]
Inheritances =[]
VariantINFO =[]
DiseaseInfo =[]
line_list =[]
VariantINFO =[]
SZAreportCategory =[]
ReviewStar = []
Compound_heterozygous = []

# Remove the duplicated features fuction
def rm_duplicates(dupfile,nodupfile):
    with open(dupfile,'r') as f:
        next(f)
        for line in f:
            line = line.replace('\n', '')  
            temp = line.split('\t') 
            iline = line 

            iClinVar_GENEINFO = temp[colNames.index('ClinVar_GENEINFO')]
            iClinVar_GENEINFO = [s.split(':')[0] for s in iClinVar_GENEINFO.split("&")]
            if len(iClinVar_GENEINFO) > 1:
                iClinVar_GENEINFO = list(set(iClinVar_GENEINFO))
                iClinVar_GENEINFO = ";".join(iClinVar_GENEINFO)
            iSZAID = temp[colNames.index('SZAID')]
            iSymbol = temp[colNames.index('SYMBOL')]
            iGene = temp[colNames.index('Genes')]
            iDiseaseInfo = temp[colNames.index('DiseaseInfo')]
            iinh = temp[colNames.index('Inheritance')]
            imark = temp[colNames.index('Mark')]
            imim = temp[colNames.index('MIM')].split(".")[0]
            iSZAreportCategory = temp[colNames.index('SZAreportCategory')]
            istar = temp[colNames.index('ClinVar_ReviewStar')]
            iREF=temp[colNames.index('REF')]
            iALT=temp[colNames.index('ALT')]
            iVARIANT_CLASS=temp[colNames.index('VARIANT_CLASS')]
            iClinvar_CLNSIG=temp[colNames.index('ClinVar_CLNSIG')].split("&")[0]
            iID = temp[colNames.index('ID')]
            iheritance = temp[colNames.index('Inheritances')]
            if iGene == 'Intergenic':
                iDisease = temp[colNames.index('ClinVar_CLNDN')]
            else:
                iDisease = temp[colNames.index('Disease')]
            if iDisease == '':
                sys.exit()
            iTargetGroup = temp[colNames.index('Target.group')]
            if iTargetGroup == '':
                iTargetGroup = 'Expanded_clinvar_mono_diseases'
            iCHROM = temp[colNames.index('#CHROM')]
            iPOS = temp[colNames.index('POS')]
            iZygosity = temp[colNames.index('Zygosity')]
            icompound_het = temp[colNames.index('Compound_heterozygous')]
            iExisting_variation = temp[colNames.index('Existing_variation')]
            iVariantINFO = iCHROM + ":" + iPOS + " " + iREF + ">" + iALT
            iConsequence = temp[colNames.index('Consequence')]
            iClinVar_CLNHGVS = temp[colNames.index('ClinVar_CLNHGVS')]
            iVariantINFO = iCHROM + ":" + iPOS + " " + iREF + ">" + iALT
            iIndex1 = [i for i, x in enumerate(TargetGroup) if  iTargetGroup == x]
            iIndex2 = [i for i, x in enumerate(Disease) if iDisease == x]
            iIndex3 = [i for i, x in enumerate(Gene) if iGene == x]
            iIndex9 = [i for i, x in enumerate(Symbol) if iSymbol == x]
            iIndex4 = [i for i, x in enumerate(SZAID) if iSZAID == x]
            iIndex5 = [i for i, x in enumerate(SZAreportCategory) if iSZAreportCategory == x]
            iIndex6 = [i for i, x in enumerate(inh) if iinh == x]
            iIndex7 = [i for i, x in enumerate(mark) if imark == x]
            iIndex8 = [i for i, x in enumerate(mim) if imim == x]
            iIndex12345678 = [x for x in iIndex1 if x in iIndex2 and x in iIndex3 and x in iIndex4 and x in iIndex5 and x in iIndex6 and x in iIndex7 and x in iIndex8 ]
            iIndex12456789 = [x for x in iIndex1 if x in iIndex2 and x in iIndex8 and x in iIndex4 and x in iIndex5 and x in iIndex6 and x in iIndex7 and x in iIndex9 ]


            if 'SZAvar' in iSZAID and (iGene in iClinVar_GENEINFO or iGene == iClinVar_GENEINFO) and ";" not in iGene:
                if len(iIndex12345678) == 0  | len(iIndex12456789) == 0:
                    append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID, iinh, imark, imim, iExisting_variation, iCHROM, iPOS,
                                    iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                    iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory, istar,icompound_het, iSymbol)

            elif 'SZAvar' in iSZAID  and ";" in iGene:
                if isinstance(iGene, list):
                    iGene = ';'.join(iGene)
                genes_set = set(iGene.split(';'))
                
                if isinstance(iClinVar_GENEINFO, list):
                    iClinVar_GENEINFO = ';'.join(iClinVar_GENEINFO)
                geneinfo_set = set(iClinVar_GENEINFO.split(';'))
                        
                if genes_set.issubset(geneinfo_set) or geneinfo_set.issubset(genes_set):
                    if len(iIndex12345678) == 0  | len(iIndex12456789) == 0:
                        append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID, iinh, imark, imim, iExisting_variation, iCHROM, iPOS,
                                    iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                    iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory, istar,icompound_het, iSymbol)

            elif 'Novel' in iSZAID:
                if len(iIndex12345678) == 0  | len(iIndex12456789) == 0 :
                    append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID,iinh, imark, imim,  iExisting_variation, iCHROM, iPOS,
                                    iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                    iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory, istar,icompound_het, iSymbol)
            
            

#write nodup file
    with open(nodupfile,'w') as f:
        colNames_new = ["final_target_group"] + colNames
        f.write("\t".join(colNames_new)+"\n")
        for j in range(0,len(Disease)):
            if TargetGroup[j] in ['Health predipositions/Disease risk','Carrier-screening','Heriditary-cancer risk syndrome','Newborn-screening'] :
                f.write('Basic (for healthy subjects),Usually used for:'+ TargetGroup[j]+'\t' + line_list[j] +"\n")
            elif TargetGroup[j] in ['Expanded_clinvar_mono_diseases','Expanded_mono_rare_diseases']:
                f.write('Extended (for potential patients),' + TargetGroup[j]+'\t' + line_list[j] + "\n")

#%%cell 9 Generate output files 
dupfile = outFile.replace('.txt','_sp_Inheritance_4.txt')
nodupfile = outFile.replace('.txt','_sp_Inheritance_4_nodup.txt')

rm_duplicates(dupfile,nodupfile) 
print("generated no dup file,", nodupfile)


#%% Cell 10 GeneBe ACMG classification
import genebe as gnb

genebe_df = pd.read_csv(nodupfile, sep="\t",header=0)
small_df = genebe_df.loc[:,["#CHROM","POS","REF","ALT"]]
small_df = small_df.rename(columns={"#CHROM":"chr","POS":"pos","REF":"ref","ALT":"alt"})
unique_small_df = small_df.drop_duplicates()
annotated_df = gnb.annotate(unique_small_df,
    genome='hg38',
    use_ensembl=False,
    use_refseq=True,
    flatten_consequences=True,
    output_format="dataframe")
annotated_df = annotated_df.rename(columns={"chr":"#CHROM","pos":"POS","ref":"REF","alt":"ALT"})

small_annotate_all = annotated_df[["#CHROM","POS","REF","ALT","gene_symbol",
                                   "acmg_score","acmg_classification",'acmg_criteria']]
small_annotate_all = small_annotate_all.rename(columns={'gene_symbol': 'Genes'})
genebe_df = pd.merge(genebe_df, small_annotate_all, how="left", on=["#CHROM","POS","REF","ALT","Genes"])


#%% Cell 11 ClinGen gene-variants classification
clingene = pd.read_csv(clingen_file, sep="\t",header=0)
genebe_df["ClinVar_ID"] = genebe_df["ClinVar_ID"].str.split("&").str[0]
clingene = clingene[['#Variation', 'ClinVar Variation Id',
      'HGNC Gene Symbol', 'Disease', 'Mondo Id',
       'Mode of Inheritance', 'Assertion', 'Applied Evidence Codes (Met)',
       'Applied Evidence Codes (Not Met)', 'Summary of interpretation']]
clingene = clingene.rename(columns={"ClinVar Variation Id": "ClinVar_ID", "HGNC Gene Symbol":"Genes",
                                    "Assertion":"ClinGen_classification","Disease":"ClinGen_disease",
                                    "#Variation":"ClinGen_variant",
                                    "Applied Evidence Codes (Met)":"ClinGen_applied_evidence_codes"})
clingen_merge_df = pd.merge(genebe_df, clingene, how="left", on=["ClinVar_ID","Genes"])



#%%cell 12 Add GenCC and ClinGen gene-disease classification
genecc_clingen_classification = pd.read_csv(geneBaseFile, sep="\t",header=0)
genecc_clingen_classification.drop(columns='GenCC_classification_clingen', inplace=True)
clingen_merge_df = pd.merge(clingen_merge_df, genecc_clingen_classification, how="left", 
                            on=['Gene.Disease.confidence.level', 'Target.group', 'Genes', 'Disease',
                                'SZAdiseaseID', 'Inheritance', 'Mark', 'MIM', ])


#%%cell 13 Add ontology  
ontology = pd.read_csv(ontology_file, sep=",",header = 0)
ontology = ontology[['Genes',"SZAdiseaseID","category"]]
clingen_merge_df = pd.merge(clingen_merge_df, ontology, on=['Genes','SZAdiseaseID'], how='left')
clingen_merge_df.to_csv(nodupfile, sep="\t",index = None, header=True) 