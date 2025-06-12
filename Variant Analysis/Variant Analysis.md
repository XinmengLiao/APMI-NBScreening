Data analysis and management of the pilot study of NBScreening were mainly conducted in R. 

#### `Major analyses.Rmd` contains the analyis of: 
  - Assigning pathogenic variants into three disease panels (childhood-onset and adulthood-onset monogenic diseases, disease carrier status) and generate the screen-positive results. 
  - Pharmacogenics (PGx) data analysis, including PGx genotypes, haplotyps basic statistics and PGx comparisons between Turkish newborns and other populations from gnomAD.
  - Discovering potential deletirous novel predicted loss-of-function variants. 
  - Figure 3A, Figure 3C, Figure 4, Figure 5 and Supplementary Figure 1A plotting.
  - Analysis for results of Supplementary Dataset Table 1 and Supplementary Dataset Table 3. 

#### `Pipeline Refinement by UKBB and rWGS.Rmd` contains the analyis of: 
  - Refining screening pipeline used in NBScreening by UK Biobank cohort.
  - Validating the Root Cause Analysis by rWGS cohort. 
    
#### `Variant statistic.Rmd` contains the analyis of: 
  - Supplementary Figure 1 and Supplementary Figure 2 plotting.
  - Identifying common and rare high-confidence loss-of-function (HC-LoF) variants from 1100 Turkish newborns, 272 Turkish healthy adults, and Turkish Variome.
  - Comparison of the common and rare HC-LoF variants among 1100 Turkish newborns, 272 Turkish healthy adults, and Turkish Variome.
  - Refining screening pipline by simulating NBScreening strategy on 272 in-house Turkish Healthy Adults. 
