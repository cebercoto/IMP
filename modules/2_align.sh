#!/bin/bash
while [[ $# > 1 ]]
do
key="$1"
case $key in
	-H|--help)
        HELP=$2
        shift
        ;;
        -l|--list)
        LIST=${2}
        shift
        ;;
        -c|--cores)
        CORES=$2
        shift
        ;;
        -b|--build)
        REFERENCE=$2
        shift
        ;;
        --pcrdup)
        PCRDUP=${2}
        shift
        ;;
        --pcrdup-alg)
        PCRDUPALG=${2}
        shift
        ;;
        --realign)
        REALIGN=${2}
        shift
        ;;
        --recalibration)
        RECALIBRATION=${2}
        shift
        ;;
        --qc)
        QC=${2}
        shift
        ;;
        --score)
        SCORE=${2}
        shift
        ;;
        -a|--align)
        ALIGN=${2}
        shift
        ;;
        -v|--variant-calling)
        VARCALL=${2}
        shift
        ;;
        --separate)
        SEPARATE=${2}
        shift
        ;;
        --delete)
        DELETE=${2}
        shift
        ;;
        --covlegend)
        COVLEGEND=${2}
        shift
        ;;
        -g|--variant-calling-algorithm)
        VARCALLALG=${2}
        shift
        ;;
        --cram)
        CRAM=${2}
        shift
        ;;
        --alt)
        ALT=${2}
        shift
        ;;
        --bed)
        BED=${2}
        shift
        ;;
        --pad)
        PAD=${2}
        shift
        ;;
        --simsam)
        SIMSAM=${2}
        shift
        ;;
        --qc-tcov)
        QCTHRESHOLDCOV=${2}
        shift
        ;;
        --qc-xaxis)
        QCXAXIS=${2}
        shift
        ;;
	--qc-alexey)
        QCALEXEY=${2}
        shift
        ;;
        --kit)
        KIT=${2}
        shift
        ;;
	--bwa-c)
	BWAC=${2}
	shift
	;;
	--bwa-r)
	BWAR=${2}
	shift
	;;
        --bampath)
        BAMPATH=${2}
        shift
        ;;
	--sample)
	SAMPLE=${2}
	shift
	;;
        --somatic)
        SOMATIC=${2}
        shift
        ;;
        --somatic-file)
        SOMATICFILE=${2}
        shift
        ;;
        --report)
        REPORT=${2}
        shift
        ;;
        --bass-maf)
        BASSMAF=${2}
        shift
        ;;
        --bass-model)
        BASSMODEL=${2}
        shift
        ;;
        --duplicates)
        DUPLICATES=${2}
        shift
        ;;
        --duplicates-file)
        DUPFILE=${2}
        shift
        ;;
        --datatype)
        DATATYPE=${2}
        shift
        ;;
esac
shift
done
i=${SAMPLE}
## Check for previous bams in the provided list so they are not aligned again
if [ -e ${i}_${REFERENCE}.bam ]; then
	echo -e `date +%y%m%d\ %H%M%S\ %s` "     Sample ${i} has already been aligned, skipping it from aligement to ${REFERENCE}" | tee -a ${LIST}.log
	vim -c "%s/${i}//e|wq" .AlignmentTemp_${LIST}/${LIST}_Temp.txt > /dev/null 2>&1
	vim -c "%s/\n\+/\r/e|wq" .AlignmentTemp_${LIST}/${LIST}_Temp.txt > /dev/null 2>&1
	echo ${i} >> .AlignmentTemp_${LIST}/${LIST}_PrevDone.txt
        rm -r .AlignmentTemp_${i}
        exit
fi
## Check that the samples fastq files exist with the correct nomencaluture
if [ ! -e `ls ${i}*fastq.gz | head -n1` ]; then
	echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Sample ${i} fastq(s) do not exist, exiting now" | tee -a ${LIST}.log
	exit
fi
## Access temp directory
cd .AlignmentTemp_${i}
## Align
echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Aligning to ${REFERENCE}" | tee -a ../${LIST}.log
if [ ${ALT} = 'true' ]; then
        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Aligning in Alt aware mode using ${REFERENCE}.fa.alt" | tee -a ../${LIST}.log
        if [ `ls ../${i}_*.fastq.gz | wc -l` -eq 2 ]; then
	        bwa mem -c ${BWAC} -r ${BWAR} -t ${CORES} -R "@RG\tID:${i}\tLB:<LIBRARY_NAME>\tSM:${i}\tPL:ILLUMINA" /data/Resources/References/${REFERENCE}/${REFERENCE}.fa ../${i}_R1_001.fastq.gz ../${i}_R2_001.fastq.gz 2> /dev/null | /data/Resources/Software/bwa.kit/k8 /data/Resources/Software/bwa.kit/bwa-postalt.js /data/Resources/References/${REFERENCE}/${REFERENCE}.fa.alt | samtools sort -O bam -l 0 -T . -o ${i}.sorted.bam > /dev/null 2>&1
        else
                bwa mem -c ${BWAC} -r ${BWAR} -t ${CORES} -R "@RG\tID:${i}\tLB:<LIBRARY_NAME>\tSM:${i}\tPL:ILLUMINA" /data/Resources/References/${REFERENCE}/${REFERENCE}.fa ../${i}.fastq.gz 2> /dev/null | /data/Resources/Software/bwa.kit/k8 /data/Resources/Software/bwa.kit/bwa-postalt.js /data/Resources/References/${REFERENCE}/${REFERENCE}.fa.alt | samtools sort -O bam -l 0 -T . -o ${i}.sorted.bam > /dev/null 2>&1
        fi
else
	bwa mem -c ${BWAC} -r ${BWAR} -t ${CORES} -R "@RG\tID:${i}\tLB:<LIBRARY_NAME>\tSM:${i}\tPL:ILLUMINA" /data/Resources/References/${REFERENCE}/${REFERENCE}.fa ../${i}_R1_001.fastq.gz ../${i}_R2_001.fastq.gz 2> /dev/null | samtools sort -O bam -l 0 -T . -o ${i}.sorted.bam > /dev/null 2>&1
fi
## Mark/Remove Duplicates
if [ ${PCRDUP} = 'remove' ]; then
        if [ ${PCRDUPALG} = 'picard' ]; then
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Marking PCR duplicates in bam file with picard" | tee -a ../${LIST}.log
                java -jar /data/Resources/Software/Javas/picard.jar MarkDuplicates I=${i}.sorted.bam O=${i}.sorted.rmdup.bam M=../${i}_${REFERENCE}.PCRDups.txt > /dev/null 2>&1
        elif [ ${PCRDUPALG} = 'samtools' ]; then
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Removing PCR duplicates from bam file with samtools" | tee -a ../${LIST}.log
                samtools rmdup ${i}.sorted.bam ${i}.sorted.rmdup.bam 2> /dev/null
                DUPREADS=$(samtools view -F 0x04 -c ${i}.sorted.rmdup.bam)
	        READS=$(samtools view -F 0x04 -c ${i}.sorted.bam)
	        DIFFDUPREADS=$(echo "${READS} - ${DUPREADS}" | bc -l )
	        PERCTDUPS=$(echo "${DIFFDUPREADS} / ${READS}" | bc -l )
	        #echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: ${DIFFDUPREADS} reads out of ${READS} were PCR duplicates (${PERCTDUPS})" | tee -a ../${LIST}.log
	        echo -e "Sample\tREAD_PAIR_DUPLICATES\tTotalReads\tPERCENT_DUPLICATION" > ${i}_${REFERENCE}.PCRDups.txt
	        echo -e "${i}\t${DIFFDUPREADS}\t${READS}\t${PERCTDUPS}" >> ${i}_${REFERENCE}.PCRDups.txt
	        cp ${i}_${REFERENCE}.PCRDups.txt ../
        else
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Algorithm not recognised, not removing/marking PCR Duplicates" | tee -a ../${LIST}.log
                mv ${i}.sorted.bam ${i}.sorted.rmdup.bam
        fi
else
	mv ${i}.sorted.bam ${i}.sorted.rmdup.bam
fi
if [ -e ${i}.sorted.bam ]; then
	rm ${i}.sorted.bam
fi
echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Indexing bam file" | tee -a ../${LIST}.log
samtools index ${i}.sorted.rmdup.bam 2> /dev/null
## Realign
if [ ${REALIGN} = 'true' ]; then
	echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Realigning" | tee -a ../${LIST}.log
	java -Xmx20g -jar /data/Resources/Software/Javas/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -o ${i}.sorted.rmdup.list -I ${i}.sorted.rmdup.bam -nt ${CORES} --allow_potentially_misencoded_quality_scores > /dev/null 2>&1
	java -Xmx20g -jar /data/Resources/Software/Javas/GenomeAnalysisTK.jar -I ${i}.sorted.rmdup.bam -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -T IndelRealigner -targetIntervals ${i}.sorted.rmdup.list -o ${i}.sorted.rmdup.realigned.bam --allow_potentially_misencoded_quality_scores --filter_bases_not_stored > /dev/null 2>&1
else
	mv ${i}.sorted.rmdup.bam ${i}.sorted.rmdup.realigned.bam
	mv ${i}.sorted.rmdup.bam.bai ${i}.sorted.rmdup.realigned.bai
fi
if [ -e ${i}.sorted.rmdup.bam ]; then
	rm ${i}.sorted.rmdup.bam
	rm ${i}.sorted.rmdup.bam.bai
fi
## Base Quality Score Recalibration
if [ ${RECALIBRATION} = 'true' ]; then
	echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Base quallity score recalibration" | tee -a ../${LIST}.log
	java -Xmx20g -jar /data/Resources/Software/Javas/GenomeAnalysisTK.jar -l INFO -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -I ${i}.sorted.rmdup.realigned.bam -T BaseRecalibrator -o ${i}.sorted.rmdup.realigned.table -knownSites /data/Resources/References/${REFERENCE}/All.vcf.gz -nct ${CORES} --allow_potentially_misencoded_quality_scores > /dev/null 2>&1 #--fix_misencoded_quality_scores
	java -Xmx20g -jar /data/Resources/Software/Javas/GenomeAnalysisTK.jar -l INFO -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -I ${i}.sorted.rmdup.realigned.bam -T PrintReads -BQSR ${i}.sorted.rmdup.realigned.table -o ${i}.sorted.rmdup.realigned.recal.bam -nct ${CORES} --allow_potentially_misencoded_quality_scores > /dev/null 2>&1
else
	mv ${i}.sorted.rmdup.realigned.bam ${i}.sorted.rmdup.realigned.recal.bam
	mv ${i}.sorted.rmdup.realigned.bai ${i}.sorted.rmdup.realigned.recal.bai
fi
if [ -e ${i}.sorted.rmdup.realigned.bam ]; then
	rm ${i}.sorted.rmdup.realigned.bam
	rm ${i}.sorted.rmdup.realigned.bai
fi
## Move final BAM to the output dir and sort out
mv ${i}.sorted.rmdup.realigned.recal.bam ../${i}_${REFERENCE}.bam
mv ${i}.sorted.rmdup.realigned.recal.bai ../${i}_${REFERENCE}.bai
#echo ${i}_${REFERENCE}.bam >> SamplesToCall.list
echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Aligment finished" | tee -a ../${LIST}.log
cd ../
rm -r .AlignmentTemp_${i}
