#!/bin/bash
#
# all submission arguments here


########################################################################################################
# PART 2: From raw basecalls to GATK-ready reads
########################################################################################################

## Recalibration
#  - GATK BQSR (Base Quality Score Recalibration) using BaseRecalibrator
#  - GATK PrintReads
#  - GATK BQSR (Base Quality Score Recalibration) using BaseRecalibrator
#  - GATK AnalyzeCovariates
## Variant Calling
#  - GATK HaplotypeCaller
#  - bgzip + tabix .g.vcf
#  - GATK ValidateVariants
## File cleanup
#  - file cleanup

########################################################################################################

# execute using: 
# bsub -J "GATKp2[1-20]" < scripts/GATK_pipeline_v2_part2.sh

source ./utils.sh

mkdir -p $root_folder/logs


####################################################
# Variables dependent on sample
####################################################


sample=$(sed "${LSB_JOBINDEX}q;d" bams_to_process.txt)
samplelog=$root_folder/logs/$sample.log

count=$(grep $sample exclude.txt | wc -l)
if [ "$count" -ge  "1" ]; then
	exit
fi



########################################################################################################
# GATK BQSR (Base Quality Score Recalibration) 
# STEP 1: Creating covariates table masking known indel, SNP sites
########################################################################################################

if [ ! -f $root_folder/logs/part_2_GATK_BQSRpre_finished_$sample.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK BaseRecalibrator Pre---" >> "$samplelog"

	time (java -Xmx12g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T BaseRecalibrator \
	-R $path_ref \
	-nct 16 \
	-I $root_folder/input_bams/$sample.reheader.bam \
	-knownSites $bundle2_8/b37/dbsnp_138.b37.vcf \
	-knownSites $bundle2_8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
	-knownSites $bundle2_8/b37/1000G_phase1.indels.b37.vcf \
	-o $root_folder/BQSR/$sample.recal_table)
	
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished GATK BaseRecalibrator Pre" >> "$samplelog" \
		&& touch $root_folder/logs/part_2_GATK_BQSRpre_finished_$sample.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), ERROR: GATK BaseRecalibrator Pre not completed. ExitValue = $exitValue" >> "$samplelog"
		exit $exitValue
	fi
else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping GATK BaseRecalibrator Pre since it was already computed***" >> "$samplelog"
fi		


########################################################################################################
# GATK BQSR (Base Quality Score Recalibration) 
# STEP 2: Writing the new BAM with recalibrated Q scores
########################################################################################################

if [ ! -f $root_folder/logs/part_2_GATK_PrintReads_finished_$sample.txt ]; then
	
	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK PrintReads---" >> "$samplelog"

	time (java -Xmx12g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T PrintReads \
	-R $path_ref \
	-nct 16 \
	-I $root_folder/input_bams/$sample.reheader.bam \
	-BQSR $root_folder/BQSR/$sample.recal_table \
	-o $root_folder/input_bams/$sample.recal.bam)

	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK PrintReads finished---" >> "$samplelog" \
		&& touch $root_folder/logs/part_2_GATK_PrintReads_finished_$sample.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), ERROR: GATK PrintReads not completed. ExitValue = $exitValue" >> "$samplelog"
		exit $exitValue
	fi

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping GATK PrintReads since it was already computed***" >> "$samplelog"
fi


####################################################
# GATK BQSR post
####################################################


if [ ! -f $root_folder/logs/part_2_GATK_BQSRpost_finished_$sample.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK BaseRecalibrator Post---" >> "$samplelog"
	
	# Count covariates in recal.BAM	to compare before/after recalibration			

	time (java -Xmx12g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T BaseRecalibrator \
	-R $path_ref \
	-nct 16 \
	-I $root_folder/input_bams/$sample.reheader.bam \
	-knownSites $bundle2_8/b37/dbsnp_138.b37.vcf \
	-knownSites $bundle2_8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
	-knownSites $bundle2_8/b37/1000G_phase1.indels.b37.vcf \
	-BQSR $root_folder/BQSR/$sample.recal_table \
	-o $root_folder/BQSR/$sample.recal_table_post)
	
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished GATK BaseRecalibrator Post" >> "$samplelog" \
		&& touch $root_folder/logs/part_2_GATK_BQSRpost_finished_$sample.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), ERROR: GATK BaseRecalibrator Post not completed. ExitValue = $exitValue" >> "$samplelog"
		exit $exitValue
	fi	
else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping GATK BaseRecalibrator Post since it was already computed***" >> "$samplelog"
fi


####################################################
# GATK AnalyzeCovariates
####################################################

if [ ! -f $root_folder/logs/part_2_GATK_AnalyzeCovariates_finished_$sample.txt ]; then
	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK AnalyzeCovariates---" >> "$samplelog"

	# Create a pdf plot to see the results of recalibration
	time (java -Xmx15g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T AnalyzeCovariates \
	-R $path_ref \
	-before $root_folder/BQSR/$sample.recal_table \
	-after $root_folder/BQSR/$sample.recal_table_post \
	-plots $root_folder/BQSR/$sample.recal_plots.pdf \
	-csv $root_folder/BQSR/$sample.recal_plots.csv \
	-l DEBUG \
	)
	
	# Catching failed java jobs ("time" doesn't change the exit value of java above)
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK AnalyzeCovariates finished---" >> "$samplelog" \
		&& touch $root_folder/logs/part_2_GATK_AnalyzeCovariates_finished_$sample.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), ERROR: GATK AnalyzeCovariates not completed. ExitValue = $exitValue" >> "$samplelog"
		exit $exitValue
	fi

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping GATK AnalyzeCovariates since it was already computed***" >> "$samplelog"
fi



####################################################
# Genotyping - HaplotypeCaller
####################################################


if [ ! -f $root_folder/logs/part_2_HaplotypeCaller_finished_$sample.txt ]; then
	echo "Starting HaplotypeCaller"
	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting HaplotypeCaller---" >> "$samplelog"

	time (java -Xmx12g -Djava.io.tmpdir=/tmp \
		-Djava.library.path=$root_folder \
		-jar $path_GATK \
		-T HaplotypeCaller \
		-R $path_ref \
		-I $root_folder/input_bams/$sample.recal.bam \
		--emitRefConfidence GVCF \
		--dbsnp $bundle2_8/b37/dbsnp_138.b37.vcf \
		-o $path_rds_vcfs/$sample.g.vcf.gz \
		-pairHMM VECTOR_LOGLESS_CACHING )
	# 		-nct 16 \

	# Catching failed java jobs ("time" doesn't change the exit value of java above)				
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished HaplotypeCaller---" >> "$samplelog" \
			&& touch $root_folder/logs/part_2_HaplotypeCaller_finished_$sample.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---HaplotypeCaller failed---" >> "$samplelog"
		exit $exitValue
	fi
	
else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping HaplotypeCaller since it was previously computed***" >> "$samplelog"
fi


####################################################
# Check the validity of the g.vcf.gz file
####################################################


if [ ! -f $root_folder/logs/part_2_GVCF_validation_finished_$sample.txt ]; then

	echo "Checking GVCF validity"
	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GVCF validation after HaplotypeCaller---" >> "$samplelog"

	time (java -Xmx12g -Djava.io.tmpdir=/tmp \
	-jar $path_GATK \
	-T ValidateVariants \
	-R $path_ref \
	-V $path_rds_vcfs/$sample.g.vcf.gz \
	--validationTypeToExclude ALL \
	) >> "$samplelog"
	# --validateGVCF \

	# NOTE that in case of invalid VCF, GATK will exit anyway

	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---GVCF validation after HaplotypeCaller COMPLETED---">> "$samplelog" \
		&& touch $root_folder/logs/part_2_GVCF_validation_finished_$sample.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---GVCF validation after HaplotypeCaller NOT COMPLETED---">> "$samplelog"
		exit $exitValue
	fi


else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping gVCF.gz validation since it was previously computed***" >> "$samplelog"
fi



####################################################
# File cleanup
####################################################


echo "$(date '+%d/%m/%y_%H:%M:%S'),---Removing last BAMs (reheader and recal)---" >> "$samplelog"
rm -f $root_folder/input_bams/$sample.reheader.*
rm -f $root_folder/input_bams/$sample.recal.*
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished removing last BAMs (reheader and recal)" >> "$samplelog"




