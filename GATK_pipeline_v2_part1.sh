#!/bin/bash
#
# all submission arguments here


########################################################################################################
# PART 1: From raw basecalls to GATK-ready reads
########################################################################################################
## Clean BAM
#  - mark duplicates using Picard tools MarkDuplicates (removed)
#  - CleanSam using Picard tools
#  - FixMateInformation using Picard tools
#  - reformat BAM header (HiSeq -> illumina)
#  - ValidateSamFile using Picard tools
#  - generate BAM index using samtools
## Metrics
#  - verify BAM id using verifyBamID
#  - get stats of BAM using Samtools flagstat
#  - wgs metrics using Picard tools CollectWgsMetrics
#  - insert size, alignment and GC bias metrics using Picard tools CollectMultipleMetrics
## File cleanup
#  - file cleanup
########################################################################################################

# execute using: 
# bsub -J "GATKp1[1-20]" < GATK_pipeline_v2_part1.sh


source ./utils.sh

mkdir -p $root_folder/logs


####################################################
# Variables dependent on sample
####################################################

# sample=$(sed "${LSB_JOBINDEX}q;d" bams_to_process_part1.txt)
sample=$(sed "${LSB_JOBINDEX}q;d" bams_to_process.txt)
samplelog=$root_folder/logs/$sample.log

count=$(grep $sample exclude.txt | wc -l)
if [ "$count" -ge  "1" ]; then
	exit
fi


####################################################
# Mark Duplicates
####################################################


# if [ ! -f $root_folder/logs/Picard_MarkDuplicates_finished_$sample.txt ]; then
# 	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting Picard MarkDuplicates---" >> "$samplelog"

# 	time (\
# 		java -Xmx2g -Djava.io.tmpdir=/tmp \
# 		-jar $path_picard MarkDuplicates \
# 		I=$original_path/$sample.bam \
# 		O=$root_folder/$sample.md.bam \
# 		M=$root_folder/logs/$sample.duplicate_metrics.txt \
# 		CREATE_INDEX=true \
# 		VALIDATION_STRINGENCY=SILENT \
# 		REMOVE_DUPLICATES=true \
# 		)


# 	exitValue=$?
# 	if [ $exitValue == 0 ]; then
# 		echo "$(date '+%d/%m/%y_%H:%M:%S'),Finished Picard MarkDuplicates" >> "$samplelog" \
# 		&& touch $root_folder/logs/Picard_MarkDuplicates_finished_$sample.txt
# 	else
# 		echo "$(date '+%d/%m/%y_%H:%M:%S'), ERROR: Picard MarkDuplicates on $sample.bam not completed. ExitValue = $exitValue" >> "$samplelog"
# 		exit $exitValue
# 	fi	

# else
# 	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping Picard MarkDuplicates since it was already computed***" >> "$samplelog"
# fi


####################################################
# Clean BAM
####################################################

if [ ! -f $root_folder/logs/part_1_Picard_CleanSam_finished_$sample.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Cleaning BAM---" >> "$samplelog"

	java -Xmx8g -Djava.io.tmpdir=/tmp \
		-jar $path_picard CleanSam \
		I=$original_path/$sample.cram \
		R=$path_ref \
		O=$root_folder/input_bams/$sample.clean.bam >> "$samplelog"

	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished cleaning BAM---" >> "$samplelog" \
		&& touch $root_folder/logs/part_1_Picard_CleanSam_finished_$sample.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), ERROR: Picard CleanSam not completed. ExitValue = $exitValue" >> "$samplelog"
		exit $exitValue
	fi

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping Picard CleanSam since it was already computed***" >> "$samplelog"
fi



####################################################
# FixMateInformation
####################################################

if [ ! -f $root_folder/logs/part_1_Picard_FixMateInformation_finished_$sample.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting FixMateInformation---" >> "$samplelog"

	time (java -Xmx8g -Djava.io.tmpdir=$tmp \
		-jar $path_picard FixMateInformation \
		I=$root_folder/input_bams/$sample.clean.bam \
		O=$root_folder/input_bams/$sample.fixed.bam \
		TMP_DIR=$tmp \
		) 

	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished Picard FixMateInformation" >> "$samplelog" \
		&& touch $root_folder/logs/part_1_Picard_FixMateInformation_finished_$sample.txt
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), ERROR: Picard FixMateInformation not completed. ExitValue = $exitValue" >> "$samplelog"
		exit $exitValue
	fi

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping FixMateInformation since it was already computed***" >> "$samplelog"
fi


####################################################
# Reformat BAM header (PL: HiSeq -> illumina)
####################################################

if [ ! -f $root_folder/logs/part_1_Reformat_header_finished_$sample.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting Reformat header---" >> "$samplelog"

	samtools view -H $root_folder/input_bams/$sample.fixed.bam | sed -e 's/PL:HiSeq/PL:illumina/g' | samtools reheader - $root_folder/input_bams/$sample.fixed.bam > $root_folder/input_bams/$sample.reheader.bam \
	&& echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished reformating header" >> "$samplelog" \
	&& touch $root_folder/logs/part_1_Reformat_header_finished_$sample.txt


else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping reformating header since it was already computed***" >> "$samplelog"
fi


####################################################
# Generate BAM index
####################################################

if [ ! -f $root_folder/input_bams/$sample.reheader.bam.bai ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Generating BAM index---" >> "$samplelog"

	time (samtools index $root_folder/input_bams/$sample.reheader.bam) \
	&& echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished generating BAM index" >> "$samplelog"

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping generating BAM index since it was already created***" >> "$samplelog"
fi


####################################################
# Validate BAM
####################################################

if [ ! -f $root_folder/logs/part_1_Picard_ValidateSamFile_finished_$sample.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Validating BAM---" >> "$samplelog"

	java -Xmx8g -Djava.io.tmpdir=/tmp \
		-jar $path_picard ValidateSamFile \
		I=$root_folder/input_bams/$sample.reheader.bam \
		MODE=SUMMARY >> "$samplelog"


	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished BAM Validation---" >> "$samplelog" \
	&& touch $root_folder/logs/part_1_Picard_ValidateSamFile_finished_$sample.txt


else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping BAM Validation since it was already computed***" >> "$samplelog"
fi


####################################################
# verifyBamID
####################################################

if [ ! -f $root_folder/logs/part_1_VerifyBamID_finished_$sample.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting VerifyBamID---" >> "$samplelog"

	time (verifyBamID \
	--vcf $vcf_REF \
	--bam $root_folder/input_bams/$sample.reheader.bam \
	--out $root_folder/logs/$sample.VerifyBamID \
	--ignoreRG \
	--verbose)

	echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished VerifyBamID" >> "$samplelog" \
	&& touch $root_folder/logs/part_1_VerifyBamID_finished_$sample.txt


else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping VerifyBamID since it was already computed***" >> "$samplelog"
fi


####################################################
# SAM FILE FLAG STATISTICS (samtools flagstat)
####################################################


if [ ! -f $root_folder/logs/part_1_Samtools_flagstats_finished_$sample.txt ]; then
	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting Samtools flagstat---" >> "$samplelog"
	samtools flagstat $root_folder/input_bams/$sample.reheader.bam >> "$samplelog" \
	&& touch $root_folder/logs/part_1_Samtools_flagstats_finished_$sample.txt \
	&& echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished samtools flagstat" >> "$samplelog"
else 
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping samtools flagstat since it was already computed***" >> "$samplelog"
fi


####################################################
# INSERT SIZE METRICS
####################################################


# if [ ! -f $root_folder/logs/part_1_Picard_CollectInsertSizeMetrics_finished_$sample.txt ]; then
# 	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting Picard CollectInsertSizeMetrics---" >> "$samplelog"

# 	time (java -Xmx8g -Djava.io.tmpdir=/tmp \
# 		-jar $path_picard CollectInsertSizeMetrics \
# 		I=$root_folder/input_bams/$sample.reheader.bam \
# 		O=$root_folder/logs/$sample.insertsize_metrics.txt \
# 		METRIC_ACCUMULATION_LEVEL=null \
# 		METRIC_ACCUMULATION_LEVEL=READ_GROUP \
# 		VALIDATION_STRINGENCY=SILENT \
# 		H=$root_folder/logs/$sample.insertsize_graph.pdf) 

# 	exitValue=$?
# 	if [ $exitValue == 0 ]; then
# 		echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished Picard CollectInsertSizeMetrics" >> "$samplelog" \
# 		&& touch $root_folder/logs/part_1_Picard_CollectInsertSizeMetrics_finished_$sample.txt
# 	else
# 		echo "$(date '+%d/%m/%y_%H:%M:%S'), ERROR: Picard CollectInsertSizeMetrics not completed. ExitValue = $exitValue" >> "$samplelog"
# 		exit $exitValue
# 	fi

# else
# 	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping Picard CollectInsertSizeMetrics since it was already computed***" >> "$samplelog"
# fi


####################################################
# CollectWgsMetrics
####################################################


if [ ! -f $root_folder/logs/part_1_Picard_CollectWgsMetrics_finished_$sample.txt ]; then
	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting Picard CollectWgsMetrics---" >> "$samplelog"

	time (java -Xmx12g -Djava.io.tmpdir=/tmp \
		-jar $path_picard CollectWgsMetrics \
		R=$path_ref \
		I=$root_folder/input_bams/$sample.reheader.bam \
		O=$root_folder/logs/$sample.wgs_metrics.txt \
		INCLUDE_BQ_HISTOGRAM=true)

	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished Picard CollectWgsMetrics" >> "$samplelog" \
		&& touch $root_folder/logs/part_1_Picard_CollectWgsMetrics_finished_$sample.txt				
	else
		echo "$(date '+%d/%m/%y_%H:%M:%S'), ERROR: Picard CollectWgsMetrics not completed. ExitValue = $exitValue" >> "$samplelog"
		exit $exitValue
	fi	

else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping Picard CollectWgsMetrics since it was already computed***" >> "$samplelog"
fi


##################################################
# Collect multiple metrics
##################################################

# Collect metrics on insert size, GC bias, alignment summary

if [ ! -f $root_folder/logs/part_1_Picard_CollectMultiMetrics_finished_$sample.txt ]; then

	echo "$(date '+%d/%m/%y_%H:%M:%S'),---Collecting multiple metrics---" >> "$samplelog"



	time (java -Xmx12g -Djava.io.tmpdir=/tmp \
		-jar $path_picard CollectMultipleMetrics \
		R=$path_ref \
		I=$root_folder/input_bams/$sample.reheader.bam \
		O=$root_folder/logs/$sample.multiple_metrics \
		PROGRAM=null \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics \
		PROGRAM=CollectGcBiasMetrics \
		METRIC_ACCUMULATION_LEVEL=null \
		METRIC_ACCUMULATION_LEVEL=READ_GROUP \
		METRIC_ACCUMULATION_LEVEL=SAMPLE \
		VALIDATION_STRINGENCY=SILENT)
			
	exitValue=$?
	if [ $exitValue == 0 ]; then
		echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished collecting multiple metrics" >> "$samplelog" \
		&& touch $root_folder/logs/part_1_Picard_CollectMultiMetrics_finished_$sample.txt
	else
		exit $exitValue
	fi	
else
	echo "$(date '+%d/%m/%y_%H:%M:%S'),***Skipping collecting multiple metrics since it was already computed***" >> "$samplelog"
fi




####################################################
# File cleanup
####################################################

echo "$(date '+%d/%m/%y_%H:%M:%S'),---Removing intermediate files---" >> "$samplelog"
rm -f $root_folder/input_bams/$sample.bam
rm -f $root_folder/input_bams/$sample.clean.bam
rm -f $root_folder/input_bams/$sample.fixed.bam
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished removing intermediate files" >> "$samplelog"













