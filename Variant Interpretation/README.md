In NBScreening, variants are called by DRAGEN. Variant annotation and interpretation were conducted by VEP (e113) for small variants and AnnotSV for structural variants. Output from VEP is converted from vcf file to txt file in Python3. 

#### `SNVs and INDELs.sh` contains the analyis of: 
  - Running VEP for SNVs and INDELs. 

#### `CNV_AnnotSV.Rmd` contains the analyis of: 
  - Running AnnotSV for CNVs. 
    
#### `vep management.py` contains the analyis of: 
  - Managing VEP outputs (vcf files) and generating tabular files. 
