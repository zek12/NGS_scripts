#!/bin/bash
#
#BSUB -J "VEP[1-25]"
# all other submission arguments here


# 1..22, 23 = X, 24 = Y, 25 = MT
v=($(seq 22))
v+=('X')
v+=('Y')
i=$((${LSB_JOBINDEX} - 1))
chr=${v[$i]}


input_file=/path_to/input.vcf.gz
output_folder="/path_to/output"

mkdir -p $output_folder/temporary_files
cd $output_folder/temporary_files

# 1) split VCF into chromosomes
bcftools view -O v -o $chr.vcf $input_file -r chr$chr

# 2) split multi-allelic into biallelic
bcftools norm -m - $chr.vcf -O v -o $chr.biallelic.vcf

# 3) remove ALT = *
bcftools view -h $chr.biallelic.vcf > header.$chr
bcftools view -H $chr.biallelic.vcf | grep -v -P '\t\*\t' > body.$chr
cat header.$chr body.$chr > $chr.biallelic.no_star.vcf
rm header.$chr body.$chr

# 4) set IDs as chr_pos_ref_alt
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' $chr.biallelic.no_star.vcf > $chr.biallelic.no_star.with_ID.vcf

# 5) add INFO/GRCh38_pos and rsID
# 5.1) First, copy current ID in INFO/GRCh38_pos
bcftools view -H $chr.biallelic.no_star.with_ID.vcf | awk '{printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $2, $4, $5, $3)}' > $chr.annotation_GRCh38.txt
bgzip $chr.annotation_GRCh38.txt -f
tabix -s 1 -b 2 -e 3 $chr.annotation_GRCh38.txt.gz

bgzip $chr.biallelic.no_star.with_ID.vcf -f
tabix -p vcf $chr.biallelic.no_star.with_ID.vcf.gz

bcftools annotate \
	--annotations $chr.annotation_GRCh38.txt.gz \
	--columns CHROM,FROM,TO,REF,ALT,GRCh38_pos \
	-h header_lines_GRCh38.txt \
	--output $chr.biallelic.no_star.with_ID.GRCh38_pos.vcf.gz --output-type z \
	$chr.biallelic.no_star.with_ID.vcf.gz

tabix -p vcf $chr.biallelic.no_star.with_ID.GRCh38_pos.vcf.gz

# 5.2) add rsID
bcftools annotate \
	--output $chr.biallelic.no_star.with_ID.GRCh38_pos.no_ids.vcf.gz \
	--output-type z \
	--remove ID \
	$chr.biallelic.no_star.with_ID.GRCh38_pos.vcf.gz
tabix -p vcf $chr.biallelic.no_star.with_ID.GRCh38_pos.no_ids.vcf.gz

bcftools annotate \
	--annotations /path_to_dbSNP_rsID/All_20180418.vcf.gz \
	--columns ID \
	--output $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.vcf \
	--output-type v \
	$chr.biallelic.no_star.with_ID.GRCh38_pos.no_ids.vcf.gz

rm $chr.biallelic.no_star.with_ID.*

# 6) flag repeated regions
java -jar /apps/gatk/3.7-0/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R /path_to_ref_genome/ref.fa \
	-V $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.vcf \
	-mask /path_to_repeat_regions/all_repeats.sorted.bed \
	-maskName 'repeat_region' \
	-o $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.repeat_flagged.vcf

# GATK VariantFiltration changes, in the VCF header, "##FORMAT=<ID=AD,Number=., ..." to "##FORMAT=<ID=AD,Number=R, ..."
# Change it back to Number=., otherwise PyVCF will complain when trying to read the VCF.
sed -i -e 's/##FORMAT=<ID=AD,Number=R/##FORMAT=<ID=AD,Number=./g' $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.repeat_flagged.vcf



# 7) do the liftover to GRCh37
# run liftOver.R to generate $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.repeat_flagged.lifted_rtracklayer.vcf


# 8) annotate using VEP
# GRCh37
/path_to_vep/ensembl-vep/vep \
	-i $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.repeat_flagged.lifted_rtracklayer.vcf \
	--no_stats \
	-o $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.repeat_flagged.lifted_rtracklayer.VEP.vcf --vcf \
	--dir_plugins /path_to_vep/ensembl-vep/loftee \
	--plugin LoF,human_ancestor_fa:/path_to_vep/ensembl-vep/loftee/human_ancestor.fa.rz,filter_position:0.05 \
	--plugin ExAC,/path_to_ExAC_mafs/ExAC.r1.sites.vep.with_chr.vcf.gz \
	--plugin CADD,/path_to_CADD_scores/whole_genome_SNVs.tsv.gz,/path_to_CADD_scores/InDels.tsv.gz \
	--plugin REVEL,/path_to_REVEL_scores/new_tabbed_revel.tsv.gz \
	--cache --force_overwrite --poly b --sift b --gene_phenotype --symbol --pick_allele --port 3337


# add --port 3337 to use Assembly GRCh37 for VEP
# add --check_existing for ClinVar
# /path_to_vep/ensembl-vep/loftee/human_ancestor.fa.rz,filter_position:0.05 is GRCh37
# REVEL scores are in GRCh37


# add CLINVAR
bgzip $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.repeat_flagged.lifted_rtracklayer.VEP.vcf
tabix -p vcf $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.repeat_flagged.lifted_rtracklayer.VEP.vcf.gz
bcftools annotate \
	--annotations /path_to_clinvar_annotations/clinvar_20181028.with_chr.vcf.gz \
	--columns "INFO/CLNSIG" \
	--output $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.repeat_flagged.lifted_rtracklayer.VEP.clinvar.vcf --output-type v \
	$chr.biallelic.no_star.GRCh38_pos.rsIDannotated.repeat_flagged.lifted_rtracklayer.VEP.vcf.gz


# 9) remove field FORMAT/PL from VCF, as in the filtering step, PyVCF returns an error when reading the VCF (cannot convert empty PL to Float)
bcftools annotate -x FORMAT/PL $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.repeat_flagged.lifted_rtracklayer.VEP.clinvar.vcf -O v -o $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.repeat_flagged.lifted_rtracklayer.VEP.clinvar.noPL.vcf


# 10) concat all chromosomes
# bcftools concat -O v -o ../annotated.vcf \
# 1.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 2.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 3.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 4.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 5.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 6.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 7.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 8.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 9.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 10.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 11.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 12.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 13.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 14.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 15.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 16.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 17.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 18.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 19.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 20.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 21.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# 22.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# X.biallelic.no_star.annotated.repeat_flagged.noPL.vcf \
# Y.biallelic.no_star.annotated.repeat_flagged.noPL.vcf

# python filter_variants.py annotated.vcf


# bsub < annotate_variants.sh


