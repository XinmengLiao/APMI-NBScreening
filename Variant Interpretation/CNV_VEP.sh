conda activate vep


zgrep -E "^#|PASS" $cnv/sample.vcf.gz | bgzip > sample_PASS.vcf.gz 


echo -e "Start running vep annotation for CNV file \ninput: ${1%%.vcf*}_PASS.vcf.gz at $(date).\n"
echo -e "vep cache: 113 \nvep version: 113 \n"

vep --cache --dir_cache vep_cache \
--cache_version 113 \
--fork 20 \
--format vcf \
--dir_plugins vep_cache/VEP_plugins/ \
-i sample_PASS.vcf.gz \
-o sample_cnv_vep.vcf.gz \
--force_overwrite \
--assembly GRCh38 \
--compress_output bgzip \
--symbol --vcf --check_existing --variant_class \
--pick --pick_order mane_select,rank \
--hgvs --refseq \
--af --af_gnomade --af_gnomadg --max_af \
--fasta vep_cache/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
--canonical \
--plugin StructuralVariantOverlap,file=vep_cache/vep_custom/nstd102.GRCh38.variant_call.vcf.gz,percentage=100,reciprocal=1,label=dbVar,cols=SVTYPE:CLNSIG:DESC:PHENO:ORIGIN \
--plugin StructuralVariantOverlap,file=vep_cache/gnomad4/gnomad.v4.1.cnv.all.vcf.gz,percentage=100,reciprocal=1,label=gnomad4e_rare_cnv,cols=SVTYPE:Genes:SF \
--plugin StructuralVariantOverlap,file=vep_cache/gnomad4/gnomad.v4.1.sv.sites.vcf.gz,percentage=100,reciprocal=1,label=gnomad4g,cols=CN_FREQ:CN_FREQ_afr:CN_FREQ_ami:CN_FREQ_asj:CN_FREQ_eas:CN_FREQ_fin:CN_FREQ_mid:CN_FREQ_nfe:CN_FREQ_rmi:CN_FREQ_sas \
--verbose

echo -e "CNV VEP finished at $(date)."


