#!/bin/bash
#
# all submission arguments here

###########################################################################################################################################################
# PART 3: Joint-genotyping
###########################################################################################################################################################
## 
# - (only when we have more than 200 samples) Combine all samples into 1 batch using CombineGVCFs (batch1.g.vcf, batch2.g.vcf, ..., batchM.g.vcf)
# - validate variants
## 
# - Joint-genotyping of all samples per chromosome (joint_chr1.vcf, joint_chr2.vcf, ..., joint_chr22.vcf, joint_chrX.vcf, joint_chrY.vcf, joint_chrMT.vcf)
# - validate variants
###########################################################################################################################################################



# execute using: 
# bsub < scripts/GATK_pipeline_v2_part3.sh


source ./utils.sh

mkdir -p $root_folder/output_gatk
mkdir -p $root_folder/logs

logfile=$root_folder/logs/log.log

# 1..22, 23 = X, 24 = Y, 25 = MT
v=($(seq 22))
v+=('X')
v+=('Y')
v+=('MT')
i=$((${LSB_JOBINDEX} - 1))
chr=${v[$i]}





####################################################
# Joint Genotyping - GenotypeGVCFs
####################################################

if [ ! -f $root_folder/logs/part_3_GenotypeGVCFs_finished_chr$chr.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GenotypeGVCFs: joint genotyping of chromosome $chr---" >> "$logfile"

	# gvcf_paths=$(for i in $(cat $root_folder/bams_to_process_part3.txt); do \
	# find $path_rds_vcfs -iname "$i".g.vcf.gz; done)
	gvcf_paths=$(ls $path_rds_vcfs/*.g.vcf.gz)

	# Create an array variable that contains location of all g.vcf input files
	gvcf_array=$(for i in $gvcf_paths; do echo "--variant $i"; done)

	# NOTE: WHEN WE HAVE MORE THAN 200 SAMPLES GVCFS, WE'LL HAVE BATCHES, SO REPLACE PREVIOUS LINES BY:
	# gvcf_array="--variant batch1.g.vcf --variant batch2.g.vcf"


	time(java -Xmx32g -Djava.io.tmpdir=/tmp \
		-jar $path_GATK \
		-T GenotypeGVCFs \
		-R $path_ref \
		--dbsnp $bundle2_8/b37/dbsnp_138.b37.vcf \
		${gvcf_array[@]} \
		-L $chr \
		-newQual \
		--disable_auto_index_creation_and_locking_when_reading_rods \
		-o $root_folder/output_gatk/joint_chr$chr.vcf )

	# added -newQual
	# if the previous command gives an error, try:
	# - removing -nt 16 \ Da error MESSAGE: Code exception (see stack trace for error itself)
	# - extracting .g.vcf.gz to .g.vcf


	# Catching failed java jobs ("time" doesn't change the exit value of java above)
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),Finished GenotypeGVCFs of chromosome $chr" >> "$logfile" \
			&& touch $root_folder/logs/part_3_GenotypeGVCFs_finished_chr$chr.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), GenotypeGVCFs of chromosome $chr NOT FINISHED" >> "$logfile"
		exit $exitValue
	fi
	
else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping GenotypeGVCFs of chromosome $chr since it was previously computed***" >> "$logfile"
fi




####################################################
# Check the validity of the vcf file
####################################################


if [ ! -f $root_folder/logs/part_3_joint_validation_finished_chr$chr.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Validating joint VCF on chromosome $chr---" >> "$logfile"

	time (java -Xmx32g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T ValidateVariants \
	-R $path_ref \
	-V $root_folder/output_gatk/joint_chr$chr.vcf \
	--validationTypeToExclude ALL \
	) >> "$logfile"
	

	# NOTE that in case of invalid VCF, GATK will exit anyway

	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), Validation of VCF chromosome $chr completed">> "$logfile" \
		&& touch $root_folder/logs/part_3_joint_validation_finished_chr$chr.txt 
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), Validation of VCF chromosome $chr NOT FINISHED">> "$logfile"
		exit $exitValue
	fi


else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping VCF validation on chromosome $chr since it was previously computed***" >> "$logfile"
fi

