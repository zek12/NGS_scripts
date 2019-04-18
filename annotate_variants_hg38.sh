module load tabix/0.2.6
module load vcftools/0.1.14
module load bcftools/1.3
module load samtools/1.5


# bsub -J "VEP[1-25]" <all_other_submission_parameters> \
# "sh annotate_variants_hg38.sh \
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
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 1/12: Extracting chromosome $chr ---"
if (( $(vcftools --gzvcf $input_file --chr $chr --stdout --recode --recode-INFO-all | bcftools view -H | wc -l) == 0 ))
then
	# try extracting including 'chr'
	if (( $(vcftools --gzvcf $input_file --chr chr$chr --stdout --recode --recode-INFO-all | bcftools view -H | wc -l) == 0 ))
	then
		echo "INFO: input VCF does not contain chromosome $chr nor chr$chr. Finishing..."
		exit
	else
		vcftools --gzvcf $input_file --chr chr$chr --out $chr --recode --recode-INFO-all # this generates 17.recode.vcf
		mv $chr.recode.vcf $chr.vcf
		# remove 'chr' notation from now on
		# awk '{gsub(/^chr/,""); print}' $chr.recode.vcf | awk '{ gsub(/^##contig=<ID=chr/,"##contig=<ID="); print }' > $chr.vcf
		# rm $chr.recode.vcf
	fi
else
	vcftools --gzvcf $input_file --chr $chr --out $chr --recode --recode-INFO-all # this generates 17.recode.vcf
	mv $chr.recode.vcf $chr.vcf
fi


# 2) sort VCF by position
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 2/12: Sorting VCF ---"
vcf-sort $chr.vcf > $chr.sorted.vcf


# 3) split multi-allelic into biallelic
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 3/12: Splitting multi-allelic into biallelic ---"
bcftools norm -m - $chr.sorted.vcf -O v -o $chr.biallelic.vcf


# 4) remove ALT = *
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 4/12: Removing ALT=* ---"
bcftools view -h $chr.biallelic.vcf > header.$chr
bcftools view -H $chr.biallelic.vcf | grep -v -P '\t\*\t' > body.$chr
cat header.$chr body.$chr > $chr.biallelic.no_star.vcf
rm header.$chr body.$chr


# 5) set IDs as CHROM_POS_REF_ALT
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 5/12: Setting ID as CHROM_POS_REF_ALT ---"
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' $chr.biallelic.no_star.vcf > $chr.biallelic.no_star.with_ID.vcf


# 6) copy current ID in INFO/GRCh38_pos
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 6/12: Adding INFO/GRCh38_pos ---"
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


# 7) add rsID
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 7/12: Adding rsIDs ---"
bcftools annotate \
	--output $chr.biallelic.no_star.with_ID.GRCh38_pos.no_ids.vcf.gz \
	--output-type z \
	--remove ID \
	$chr.biallelic.no_star.with_ID.GRCh38_pos.vcf.gz
tabix -p vcf $chr.biallelic.no_star.with_ID.GRCh38_pos.no_ids.vcf.gz

bcftools annotate \
	--annotations /path_to_reference_files/dbSNP_rsids/GRCh38/All_20180418.vcf.gz \
	--columns ID \
	--output $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.vcf \
	--output-type v \
	$chr.biallelic.no_star.with_ID.GRCh38_pos.no_ids.vcf.gz

rm $chr.biallelic.no_star.with_ID.*


# 8) flag repeated regions
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 8/12: Flagging repeat regions ---"

# optional: correct contig lengths
# useful links about sed:
# https://unix.stackexchange.com/questions/32908/how-to-insert-the-content-of-a-file-into-another-file-before-a-pattern-marker
# https://askubuntu.com/questions/76808/how-do-i-use-variables-in-a-sed-command
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- Update contigs in VCF using $ref_genome ---"
awk '{ printf "##contig=<ID=%s,length=%d>\n", $1, $2 }' $ref_genome.fai > $chr.contigs
sed '/##contig=/d' $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.vcf > $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.vcf
sed -n -i -e "/^#CHROM/r $chr.contigs" -e 1x -e '2,${x;p}' -e '${x;p}' $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.vcf

echo "$(date '+%d/%m/%y_%H:%M:%S'),--- Running GATK VariantFiltration ---"
java -jar GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R $ref_genome \
	-V $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.vcf \
	-mask /path_to_reference_files/repeat_regions/GRCh38/all_repeats.sorted.bed \
	-maskName 'repeat_region' \
	-o $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.vcf

# GATK VariantFiltration changes, in the VCF header, "##FORMAT=<ID=AD,Number=., ..." to "##FORMAT=<ID=AD,Number=R, ..."
# Change it back to Number=., otherwise PyVCF will complain when trying to read the VCF.
sed -i -e 's/##FORMAT=<ID=AD,Number=R/##FORMAT=<ID=AD,Number=./g' $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.vcf



# 9) add ClinVar annotation in INFO/CLNSIG
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 9/12: Adding ClinVar INFO/CLNSIG ---"
bgzip $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.vcf -f
tabix -p vcf $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.vcf.gz -f
bcftools annotate \
	--annotations /path_to_reference_files/ClinVar/GRCh38/clinvar_20181028.with_chr.vcf.gz \
	--columns "INFO/CLNSIG" \
	--output $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.vcf --output-type v \
	$chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.vcf.gz


# 10) annotate using VEP
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 10/12: Waiting... Do the the liftover to hg19 first ---"
# run liftOver.R to generate $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.vcf

exit

# 11) annotate using VEP
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 11/12: Annotating with VEP ---"

# in GRCh38
# PROBLEM: doesn't have REVEL scores in GRCh38

# /path_to_vep/VEPv96/ensembl-vep/vep \
# 	-i $chr.biallelic.no_star.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.no_VEP.vcf \
# 	-o $chr.biallelic.no_star.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.VEP.vcf --vcf \
# 	--no_stats --force_overwrite --cache --pick_allele --transcript_version

# /path_to_vep/VEPv96/ensembl-vep/vep \

# /path_to_vep/ensembl-vep/vep \
# 	-i $chr.biallelic.no_star.repeat_flagged.clinvar.no_VEP.vcf \
# 	--no_stats \
# 	-o $chr.biallelic.no_star.repeat_flagged.clinvar.VEP.vcf --vcf \
# 	--dir_plugins /path_to_vep/ensembl-vep/loftee \
# 	--plugin LoF,human_ancestor_fa:/path_to_vep/ensembl-vep/loftee/human_ancestor.fa.rz,filter_position:0.05 \
# 	--plugin ExAC,/path_to_reference_files/ExAC_MAFs/GRCh38/legacy-exacv1_downloads-liftover_grch38-release1-ExAC.r1.sites.liftover.b38.vcf.gz \
# 	-custom /path_to_reference_files/gnomAD_MAFs/GRCh38/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz,gnomAD_exomes,vcf,exact,0,AF_NFE \
# 	--plugin CADD,/path_to_reference_files/CADD_scores/GRCh38/whole_genome_SNVs.tsv.gz,/path_to_reference_files/CADD_scores/GRCh38/InDels.tsv.gz \
# 	--cache --force_overwrite --poly b --sift b --gene_phenotype --symbol --pick_allele --port 3337


# remove 'chr' notation from now on
awk '{gsub(/^chr/,""); print}' $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.vcf | awk '{ gsub(/^##contig=<ID=chr/,"##contig=<ID="); print }' > \
$chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.no_chr.vcf

# remove INFO/CSQ
if (( $(bcftools view -h $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.no_chr.vcf | grep "##INFO=<ID=CSQ" | wc -l) == 0 ))
then
	cp $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.no_chr.vcf $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.no_chr.no_VEP.vcf
else	
	bcftools annotate -x INFO/CSQ $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.no_chr.vcf -O v -o $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.no_chr.no_VEP.vcf
fi

# GRCh37
/path_to_vep/ensembl-vep/vep \
	-i $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.no_chr.no_VEP.vcf \
	--no_stats \
	-o $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.VEP.vcf --vcf \
	--dir_plugins /path_to_vep/ensembl-vep/loftee \
	--plugin LoF,human_ancestor_fa:/path_to_vep/ensembl-vep/loftee/human_ancestor.fa.rz,filter_position:0.05 \
	--plugin ExAC,/path_to_reference_files/ExAC_MAFs/GRCh37/ExAC.r1.sites.vep.vcf.gz \
	-custom /path_to_reference_files/gnomAD_MAFs/GRCh37/gnomad.exomes.r2.1.sites.vcf.gz,gnomAD_exomes,vcf,exact,0,AF_nfe \
	--plugin CADD,/path_to_reference_files/CADD_scores/GRCh37/whole_genome_SNVs.tsv.gz,/path_to_reference_files/CADD_scores/GRCh37/InDels.tsv.gz \
	--plugin REVEL,/path_to_reference_files/REVEL_scores/new_tabbed_revel.tsv.gz \
	--cache --force_overwrite --poly b --sift b --gene_phenotype --symbol --pick_allele --port 3337

rm $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.no_chr.*

# more info: https://www.ensembl.org/info/docs/tools/vep/script/vep_example.html#gnomad, section "gnomAD and ExAC"
# --af_gnomad adds gnomAD exomes MAFs v2.0.1, in INFO/CSQ (gnomAD_NFE_AF)
# -custom adds gnomAD exomes MAFs v2.1, in INFO/CSQ, as the last value (gnomAD_exomes_AF_nfe)

# add --port 3337 to use Assembly GRCh37 for VEP
# add --check_existing for ClinVar
# /path_to_vep/ensembl-vep/loftee/human_ancestor.fa.rz,filter_position:0.05 is GRCh37
# REVEL scores are in GRCh37


# 12) remove field FORMAT/PL from VCF, as in the filtering step, PyVCF returns an error when reading the VCF (cannot convert empty PL to Float)
echo "$(date '+%d/%m/%y_%H:%M:%S'),--- STEP 12/12: Removing FORMAT/PL ---"
bcftools annotate -x FORMAT/PL $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.VEP.vcf -O v -o $chr.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.VEP.noPL.vcf


echo "$(date '+%d/%m/%y_%H:%M:%S'),--- DONE ---"

# 9) concat all chromosomes
# bcftools concat -O v -o ../final.vcf \
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
# bcftools concat -O v -o ../final.vcf \
# *.noPL.vcf
# cd ..
# bgzip final.vcf
# tabix -p vcf final.vcf.gz
# java -jar picard.jar SortVcf I=final.vcf.gz O=final.sorted.vcf 




