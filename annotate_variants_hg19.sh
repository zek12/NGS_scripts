module load tabix/0.2.6
module load vcftools/0.1.14
module load bcftools/1.3
module load samtools/1.5


# bsub -J "VEP[1-25]" <all_other_submission_parameters> \
# "sh annotate_variants_hg19.sh \
# -i input.vcf \
# -o output_dir/ \
# -R ref.fa"


# 1..22, 23 = X, 24 = Y, 25 = MT
v=($(seq 22))
v+=('X')
v+=('Y')
v+=('MT')
i=$((${LSB_JOBINDEX} - 1))
chr=${v[$i]}



while getopts i:o:g:R: option 
do 
	case "${option}" 
	in 
	i) input_file=${OPTARG};;
	o) output_folder=${OPTARG};;
	R) ref_genome=${OPTARG};;
	esac 
done

echo ""
echo "////////////////////////////////////////////////////"
echo "Input VCF:" $input_file
echo "Output folder:" $output_folder
echo "Reference Genome:" $ref_genome
echo "////////////////////////////////////////////////////"
echo ""


if [ ! -f $input_file ]
then
	echo "ERROR: Input VCF does not exist"
	exit
fi

if [ ! -d $output_folder ]
then
	echo "ERROR: Output folder does not exist"
	exit
fi

if [ ! -f $ref_genome ]
then
	echo "ERROR: Reference Genome does not exist"
	exit
fi

mkdir -p $output_folder/temporary_files
cd $output_folder/temporary_files



# 1) check if input VCF contains chromosome: if so, extract specific chromosome; otherwise, exit.
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 1/7: Extracting chromosome $chr ---"
if (( $(vcftools --gzvcf $input_file --chr $chr --stdout --recode --recode-INFO-all | bcftools view -H | wc -l) == 0 ))
then
	# try extracting including 'chr'
	if (( $(vcftools --gzvcf $input_file --chr chr$chr --stdout --recode --recode-INFO-all | bcftools view -H | wc -l) == 0 ))
	then
		echo "INFO: input VCF does not contain chromosome $chr nor chr$chr. Finishing..."
		exit
	else
		vcftools --gzvcf $input_file --chr chr$chr --out $chr --recode --recode-INFO-all # this generates 17.recode.vcf
		# remove 'chr' notation from now on
		awk '{gsub(/^chr/,""); print}' $chr.recode.vcf | awk '{ gsub(/^##contig=<ID=chr/,"##contig=<ID="); print }' > $chr.vcf
		rm $chr.recode.vcf
	fi
else
	vcftools --gzvcf $input_file --chr $chr --out $chr --recode --recode-INFO-all # this generates 17.recode.vcf
	mv $chr.recode.vcf $chr.vcf
fi


# 2) sort VCF by position
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 2/7: Sorting VCF ---"
vcf-sort $chr.vcf > $chr.sorted.vcf


# 3) split multi-allelic into biallelic
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 3/7: Splitting multi-allelic into biallelic ---"
bcftools norm -m - $chr.sorted.vcf -O v -o $chr.biallelic.vcf

# 4) remove ALT = *
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 4/7: Removing ALT=* ---"
bcftools view -h $chr.biallelic.vcf > header.$chr
bcftools view -H $chr.biallelic.vcf | grep -v -P '\t\*\t' > body.$chr
cat header.$chr body.$chr > $chr.biallelic.no_star.vcf
rm header.$chr body.$chr



# 5) flag repeated regions
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 5/7: Flagging repeat regions ---"

java -jar /apps/gatk/3.7-0/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R $ref_genome \
	-V $chr.biallelic.no_star.vcf \
	-mask /path_to_repeat_regions/all_repeats.sorted.bed \
	-maskName 'repeat_region' \
	-o $chr.biallelic.no_star.repeat_flagged.vcf

# GATK VariantFiltration changes, in the VCF header, "##FORMAT=<ID=AD,Number=., ..." to "##FORMAT=<ID=AD,Number=R, ..."
# Change it back to Number=., otherwise PyVCF will complain when trying to read the VCF.
sed -i -e 's/##FORMAT=<ID=AD,Number=R/##FORMAT=<ID=AD,Number=./g' $chr.biallelic.no_star.repeat_flagged.vcf



# 6) add ClinVar annotation in INFO/CLNSIG
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 6/7: Adding ClinVar INFO/CLNSIG ---"
bgzip $chr.biallelic.no_star.repeat_flagged.vcf -f
tabix -p vcf $chr.biallelic.no_star.repeat_flagged.vcf.gz -f
bcftools annotate \
	--annotations /path_to_clinvar_annotations/clinvar_20181028.vcf.gz \
	--columns "INFO/CLNSIG" \
	--output $chr.biallelic.no_star.repeat_flagged.clinvar.vcf --output-type v \
	$chr.biallelic.no_star.repeat_flagged.vcf.gz

# 7) annotate using VEP
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 7/7: Annotating with VEP ---"
# remove INFO/CSQ
bcftools annotate -x INFO/CSQ $chr.biallelic.no_star.repeat_flagged.clinvar.vcf -O v -o $chr.biallelic.no_star.repeat_flagged.clinvar.no_vep.vcf

# GRCh37
/path_to_vep/ensembl-vep/vep \
	-i $chr.biallelic.no_star.repeat_flagged.clinvar.no_vep.vcf \
	--no_stats \
	-o $chr.biallelic.no_star.repeat_flagged.clinvar.VEP.vcf --vcf \
	--dir_plugins /path_to_vep/ensembl-vep/loftee \
	--plugin LoF,human_ancestor_fa:/path_to_vep/ensembl-vep/loftee/human_ancestor.fa.rz,filter_position:0.05 \
	--plugin ExAC,/path_to_ExAC_mafs/ExAC.r1.sites.vep.vcf.gz \
	-custom /path_to_gnomAD_mafs/gnomad.exomes.r2.1.sites.vcf.gz,gnomAD_exomes,vcf,exact,0,AF_nfe \
	--plugin CADD,/path_to_CADD_scores/whole_genome_SNVs.tsv.gz,/path_to_CADD_scores/InDels.tsv.gz \
	--plugin REVEL,/path_to_REVEL_scores/new_tabbed_revel.tsv.gz \
	--cache --force_overwrite --poly b --sift b --gene_phenotype --symbol --pick_allele --port 3337
	

# more info: https://www.ensembl.org/info/docs/tools/vep/script/vep_example.html#gnomad, section "gnomAD and ExAC"
# --af_gnomad adds gnomAD exomes MAFs v2.0.1, in INFO/CSQ (gnomAD_NFE_AF)
# -custom adds gnomAD exomes MAFs v2.1, in INFO/CSQ, as the last value (gnomAD_exomes_AF_nfe)


# add --port 3337 to use Assembly GRCh37 for VEP
# add --check_existing for ClinVar
# /path_to_vep/ensembl-vep/loftee/human_ancestor.fa.rz,filter_position:0.05 is GRCh37
# REVEL scores are in GRCh37


# 8) remove field FORMAT/PL from VCF, as in the filtering step, PyVCF returns an error when reading the VCF (cannot convert empty PL to Float)
bcftools annotate -x FORMAT/PL $chr.biallelic.no_star.repeat_flagged.clinvar.VEP.vcf -O v -o $chr.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf


echo "$(date '+%d/%m/%y_%H:%M:%S'),--- DONE ---"

# 9) concat all chromosomes
# bcftools concat -O v -o ../annotated.vcf \
# 1.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 2.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 3.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 4.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 5.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 6.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 7.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 8.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 9.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 10.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 11.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 12.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 13.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 14.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 15.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 16.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 17.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 18.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 19.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 20.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 21.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# 22.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# X.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# Y.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf \
# MT.biallelic.no_star.repeat_flagged.clinvar.VEP.noPL.vcf

# or:
# bcftools concat -O v -o ../annotated.vcf *.noPL.vcf
# cd ..
# bgzip annotated.vcf
# tabix -p vcf annotated.vcf.gz
# sort: java -jar /apps/picard-tools/2.0.1/picard.jar SortVcf I=annotated.vcf.gz O=annotated.sorted.vcf
# sort: vcf-sort annotated.vcf.gz > annotated.sorted.vcf

# python filter_variants.py annotated.vcf

