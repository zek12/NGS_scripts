#!/bin/bash
#
# all submission arguments here

########################################################################################################
# PART 4: Joint genotyping of all samples, Variant Filtration
########################################################################################################
## CatVariants (merge all joint_<chr>.vcf into joint.vcf)
## Variant Recalibration
# - VariantRecalibrator (SNPs)
# - ApplyRecalibration (SNPs)
# - validate variants (SNPs)
# - VariantRecalibrator (indels)
# - ApplyRecalibration (indels)
# - validate variants (indels)
## Refinement
# - CalculateGenotypePosteriors
# - VariantFiltration (skipped)
# - validate variants
## clean VCF
# - LeftAlignAndTrimVariants
# - validate variants
## File cleanup
#  - file cleanup
########################################################################################################



# execute using:
# bsub < scripts/GATK_pipeline_v2_part4.sh


source ./utils.sh

mkdir -p $root_folder/output_gatk
mkdir -p $root_folder/logs


# num_samples=$( wc -l $root_folder/bams_to_process_part3.txt )
num_samples=$( ls $path_rds_vcfs/*.g.vcf.gz | wc -l )
# num_samples=$( echo $num_samples | cut -d' ' -f1  )
# output_name="joint_"$num_samples"samples"
output_name="227uk.germline_calls.noQC"





####################################################
# GATK CatVariants
####################################################

if [ ! -f $root_folder/logs/part_4_CatVariants_finished.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK CatVariants: merging all joint_<chr>.vcf into joint.vcf---" >> "$logfile"

	gvcf_paths=$( find output_gatk/ -iname joint_chr*.vcf )

	# Create an array variable that contains location of all joint_<chr>.vcf
	vcf_array=$(for i in $gvcf_paths; do echo "--variant $i"; done)


	time(java -Xmx16g -Djava.io.tmpdir=/tmp \
		-cp $path_GATK \
		org.broadinstitute.gatk.tools.CatVariants \
		-R $path_ref \
		${vcf_array[@]} \
		-out $root_folder/output_gatk/$output_name.vcf )


	# Catching failed java jobs ("time" doesn't change the exit value of java above)
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished GATK CatVariants" >> "$logfile" \
			&& touch $root_folder/logs/part_4_CatVariants_finished.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), GATK CatVariants NOT FINISHED" >> "$logfile"
		exit $exitValue
	fi
	
else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping GATK CatVariants since it was previously computed***" >> "$logfile"
fi



####################################################
# Check the validity of the vcf file
####################################################


if [ ! -f $root_folder/logs/part_4_VCF_validation_after_CatVariants_finished.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Validating VCF after CatVariants---" >> "$logfile"

	time (java -Xmx15g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T ValidateVariants \
	-R $path_ref \
	-V $root_folder/output_gatk/$output_name.vcf \
	--validationTypeToExclude ALL \
	) >> "$logfile"
	

	# NOTE that in case of invalid VCF, GATK will exit anyway

	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), VCF validation after CatVariants COMPLETED">> "$logfile" \
		&& touch $root_folder/logs/part_4_VCF_validation_after_CatVariants_finished.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), VCF validation after CatVariants NOT COMPLETED">> "$logfile"
		exit $exitValue
	fi


else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping VCF validation after CatVariants since it was previously computed***" >> "$logfile"
fi



############################################################################
# SNP Filtering (VQSR) - Variant Recalibrator
############################################################################


if [ ! -f $root_folder/logs/part_4_GATK_SNP_VQSR_finished.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK SNP VQSR---" >> "$logfile"

	# Creating the machine learning SNP calibration model

	time(java -Xmx15g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T VariantRecalibrator \
	-nt 16 \
	-R $path_ref \
	-input $root_folder/output_gatk/$output_name.vcf \
	-recalFile  $root_folder/output_gatk/$output_name.SNP.recal \
	-tranchesFile $root_folder/output_gatk/$output_name.SNP.tranches \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $bundle2_8/b37/hapmap_3.3.b37.vcf \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 $bundle2_8/b37/1000G_omni2.5.b37.vcf \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $bundle2_8/b37/1000G_phase1.snps.high_confidence.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $bundle2_8/b37/dbsnp_138.b37.vcf \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
	-rscriptFile $root_folder/output_gatk/$output_name.recalibrate_SNP_plots.R \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
	-mode SNP)

	# note:
	# added -an DP
	# added line -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \ 
	# eventually, we should have:
	# -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
	# see http://gatkforums.broadinstitute.org/gatk/discussion/1259/which-training-sets-arguments-should-i-use-for-running-vqsr

	
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished SNP VQSR VariantRecalibrator---" >> "$logfile"	
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK SNP VQSR VariantRecalibrator NOT FINISHED---" >> "$logfile" \
		exit $exitValue
	fi

	# Applying the SNP model

	time (java -Xmx15g -Djava.io.tmpdir=/tmp \
		-jar $path_GATK \
		-T ApplyRecalibration \
		-R $path_ref \
		-input $root_folder/output_gatk/$output_name.vcf \
		-tranchesFile $root_folder/output_gatk/$output_name.SNP.tranches \
		-recalFile $root_folder/output_gatk/$output_name.SNP.recal \
		-o $root_folder/output_gatk/$output_name.SNP.recal.vcf \
		--ts_filter_level 99.9 \
		-mode SNP)

	# --excludeFiltered \

	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished GATK SNP VQSR ApplyRecalibration---" >> "$logfile" \
		&& touch $root_folder/logs/part_4_GATK_SNP_VQSR_finished.txt		
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK SNP VQSR ApplyRecalibration NOT FINISHED---" >> "$logfile" \
		exit $exitValue
	fi

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping GATK SNP VQSR since it was previously computed***" >> "$logfile"
fi		


####################################################
# Check the validity of the vcf file
####################################################


if [ ! -f $root_folder/logs/part_4_GATK_SNP_VQSR_VCF_validation_finished.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Checking SNP_VQSR vcf validity---" >> "$logfile"

	time (java -Xmx15g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T ValidateVariants \
	-R $path_ref \
	-V $root_folder/output_gatk/$output_name.SNP.recal.vcf \
	--validationTypeToExclude ALL \
	) >> "$logfile"

	# NOTE that in case of invalid VCF, GATK will exit anyway
	
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), SNP_VQSR vcf validation completed">> "$logfile" \
		&& touch $root_folder/logs/part_4_GATK_SNP_VQSR_VCF_validation_finished.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), SNP_VQSR vcf validation NOT FINISHED">> "$logfile" \
		exit $exitValue
	fi


else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping SNP_VQSR vcf validation since it was previously computed***" >> "$logfile"
fi




############################################################################
# INDEL Filtering (VQSR) - Variant Recalibrator
############################################################################

if [ ! -f $root_folder/logs/part_4_GATK_indel_VQSR_finished.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK indel VQSR---" >> "$logfile"

	# Creating the machine learned INDEL calibration model
	# NOTE:
	# changed -resource:mills,known=true...
	# to -resource:mills,known=false...
	
	time(java -Xmx15g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T VariantRecalibrator \
	-nt 16 \
	-R $path_ref \
	-input $root_folder/output_gatk/$output_name.SNP.recal.vcf \
	-recalFile $root_folder/output_gatk/$output_name.indel.recal \
	-tranchesFile $root_folder/output_gatk/$output_name.indel.tranches \
	-resource:mills,known=false,training=true,truth=true,prior=12.0 $bundle2_8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $bundle2_8/b37/dbsnp_138.b37.vcf \
	-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
	--maxGaussians 4 \
	-rscriptFile $root_folder/output_gatk/$output_name.recalibrate_INDEL_plots.R \
	-mode INDEL)


	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished indel VQSR VariantRecalibrator---" >> "$logfile"	
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK indel VQSR VariantRecalibrator NOT FINISHED---" >> "$logfile"
		exit $exitValue
	fi


	# Applying the INDEL model

	time (java -Xmx15g -Djava.io.tmpdir=/tmp \
		-jar $path_GATK \
		-T ApplyRecalibration \
		-R $path_ref \
		-input $root_folder/output_gatk/$output_name.SNP.recal.vcf \
		-tranchesFile $root_folder/output_gatk/$output_name.indel.tranches \
		-recalFile $root_folder/output_gatk/$output_name.indel.recal \
		-o $root_folder/output_gatk/$output_name.SNP.indel.recal.vcf \
		--ts_filter_level 99.9 \
		-mode INDEL)

	# --excludeFiltered \


	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished GATK indel VQSR ApplyRecalibration---" >> "$logfile" \
		&& touch $root_folder/logs/part_4_GATK_indel_VQSR_finished.txt		
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK indel VQSR ApplyRecalibration NOT FINISHED---" >> "$logfile"	
		exit $exitValue
	fi

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping GATK indel VQSR since it was previously computed***" >> "$logfile"
fi	


####################################################
# Check the validity of the vcf file
####################################################

if [ ! -f $root_folder/logs/part_4_indel_VQSR_VCF_validation_finished.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Checking indel_VQSR vcf validity---" >> "$logfile"

	time (java -Xmx15g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T ValidateVariants \
	-R $path_ref \
	-V $root_folder/output_gatk/$output_name.SNP.indel.recal.vcf \
	--validationTypeToExclude ALL \
	) >> "$logfile"
		
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK indel_VQSR vcf validation completed---">> "$logfile" \
		&& touch $root_folder/logs/part_4_indel_VQSR_VCF_validation_finished.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---Validation of $root_folder/output_gatk/$output_name.SNP.indel.recal.vcf NOT FINISHED---">> "$logfile" \
		exit $exitValue
	fi


else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping indel_VQSR vcf validation since it was previously computed***" >> "$logfile"
fi






####################################################
# Genotype Refinement
####################################################


if [ ! -f $root_folder/logs/part_4_GATK_CalculateGenotypePosteriors_finished.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK CalculateGenotypePosteriors---" >> "$logfile"

	# Derive posterior probabilities of genotypes (using 1000G phase 3 data)

	time (java -Xmx15g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-R $path_ref \
	-T CalculateGenotypePosteriors \
	--supporting $bundle2_8/b37/1000G_phase3_v4_20130502.sites.vcf.gz \
	-V $root_folder/output_gatk/$output_name.SNP.indel.recal.vcf \
	-o $root_folder/output_gatk/$output_name.SNP.indel.recal.postCGP.vcf)
	
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished GATK CalculateGenotypePosteriors---" >> "$logfile"
		touch $root_folder/logs/part_4_GATK_CalculateGenotypePosteriors_finished.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK CalculateGenotypePosteriors NOT FINISHED---" >> "$logfile"
		exit $exitValue
	fi

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping GATK CalculateGenotypePosteriors since it was previously computed***" >> "$logfile"
fi



####################################################
# VariantFiltration
####################################################

# if [ ! -f $root_folder/logs/part_4_GATK_VariantFiltration_finished.txt ]; then

# 	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK VariantFiltration---" >> "$logfile"

# 	# set FILTER with GQ < 20 or DP < 8 to "low_GQ_or_DP", and set their GT to no call.

# 	time (java -Xmx15g -Djava.io.tmpdir=/tmp \
# 	-jar $path_GATK \
# 	-T VariantFiltration \
# 	-R $path_ref \
# 	-V $root_folder/output_gatk/$output_name.SNP.indel.recal.postCGP.vcf \
# 	--genotypeFilterExpression "GQ < 20.0 || DP < 8" --genotypeFilterName "low_GQ_or_DP" \
# 	--setFilteredGtToNocall \
# 	-o $root_folder/output_gatk/$output_name.SNP.indel.recal.postCGP.filtered.vcf)

# 	exitValue=$?
# 	if [ $exitValue == 0 ]; then
# 		echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished GATK VariantFiltration---" >> "$logfile" \
# 		&& touch $root_folder/logs/part_4_GATK_VariantFiltration_finished.txt
# 	else
# 		echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK VariantFiltration NOT FINISHED---" >> "$logfile"
# 		exit $exitValue
# 	fi

# else
# 	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping GATK VariantFiltration since it was previously computed***" >> "$logfile"
# fi


####################################################
# Check the validity of the vcf file
####################################################

if [ ! -f $root_folder/logs/part_4_postCGP_VCF_validation_finished.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Checking postCGP vcf validity---" >> "$logfile"

	time (java -Xmx15g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T ValidateVariants \
	-R $path_ref \
	-V $root_folder/output_gatk/$output_name.SNP.indel.recal.postCGP.vcf \
	--validationTypeToExclude ALL \
	) >> "$logfile"
	
	
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), postCGP vcf validation completed">> "$logfile" \
		&& touch $root_folder/logs/part_4_postCGP_VCF_validation_finished.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), postCGP vcf validation NOT COMPLETED">> "$logfile"
		exit $exitValue
	fi

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping postCGP vcf validation since it was previously computed***" >> "$logfile"
fi




####################################################
# LEFT ALIGN, TRIM, SPLIT VARIANTS
####################################################

# Left align and trim the variants, split multiple into biallelic (pre-processing for annotation)

if [ ! -f $root_folder/logs/part_4_leftalign_finished.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Left aligning, trimming and splitting variants---" >> "$logfile"

	time (java -Xmx15g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T LeftAlignAndTrimVariants \
	-R $path_ref \
	-V $root_folder/output_gatk/$output_name.SNP.indel.recal.postCGP.vcf \
	-o $root_folder/output_gatk/$output_name.SNP.indel.recal.postCGP.aln_trim_split.vcf \
	--splitMultiallelics)
	
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), Left aligning, trimming and splitting variants completed">> "$logfile" \
		&& touch $root_folder/logs/part_4_leftalign_finished.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), Left aligning, trimming and splitting variants NOT COMPLETED">> "$logfile"
		exit $exitValue
	fi

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping Left aligning, trimming and splitting variants since it was previously computed***" >> "$logfile"
fi


####################################################
# Check the validity of the vcf file
####################################################

if [ ! -f $root_folder/logs/part_4_postTrim_VCF_validation_finished.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Checking postTrim vcf validity---" >> "$logfile"

	time (java -Xmx15g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T ValidateVariants \
	-R $path_ref \
	-V $root_folder/output_gatk/$output_name.SNP.indel.recal.postCGP.aln_trim_split.vcf \
	--validationTypeToExclude ALL \
	) >> "$logfile"
	
	# NOTE that in case of invalid VCF, GATK will exit anyway
	
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), postTrim vcf validation completed">> "$logfile" \
		&& touch $root_folder/logs/part_4_postTrim_VCF_validation_finished.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), postTrim vcf validation NOT COMPLETED">> "$logfile"
		exit $exitValue
	fi

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping postTrim vcf validation since it was previously computed***" >> "$logfile"
fi




