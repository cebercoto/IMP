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
        --varmapq)
	VARMAPQ=${2}
	shift
	;;
        --bampath)
        BAMPATH=${2}
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
echo -e `date +%y%m%d\ %H%M%S\ %s` "     Performing QC on all processed samples" | tee -a ${LIST}.log
CORES2=`echo "${CORES} * ${SIMSAM}" | bc -l`
## Create QC Folder for each sample and Collect QC from previous steps: PCRdups (if PCRDUP=remove) from Alignment
for i in `cat ${LIST}.txt`; do
        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     ${i}: Preparing" | tee -a ${LIST}.log
        mkdir ${i}_${REFERENCE}_QC
        if [ -e ${i}_${REFERENCE}.PCRDups.txt ]; then
                mv ${i}_${REFERENCE}.PCRDups.txt ${i}_${REFERENCE}_QC/
        fi
done
wait
## fastqc reports generation
for i in `cat ${LIST}.txt`; do
        ONGOING=`ps aux | grep FastQC | wc -l`
        while [ ${ONGOING} -gt ${CORES2} ]; do
                ONGOING=`ps aux | grep FastQC | wc -l`
                echo -ne "\r${ONGOING} fastqc processes going on ATM, waiting for some to finish"
                sleep 1
        done
        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     ${i}: Generating FastQC report" | tee -a ${LIST}.log
        fastqc ${i}_${REFERENCE}.bam > /dev/null 2>&1 &
done
wait
for i in `cat ${LIST}.txt`; do
	unzip -qq ${i}_${REFERENCE}_fastqc.zip
	cp ${i}_${REFERENCE}_fastqc/summary.txt ${i}_${REFERENCE}_fastqc.txt
	vim -c "%s/ /_/g|%s/\(\S\+\)\s\+\(\S\+\).\+/\2\t\1/|1s/\(.\+\)/Sample\t${i}/|wq" ${i}_${REFERENCE}_fastqc.txt > /dev/null 2>&1
        awk '
        {
            for (i=1; i<=NF; i++)  {
                a[NR,i] = $i
            }
        }
        NF>p { p = NF }
        END {
            for(j=1; j<=p; j++) {
                str=a[1,j]
                for(i=2; i<=NR; i++){
                    str=str" "a[i,j];
                }
                print str
            }
        }' ${i}_${REFERENCE}_fastqc.txt | sed -e "s/\s\+/\t/" > ${i}_${REFERENCE}_fastqc.txt2
        mv ${i}_${REFERENCE}_fastqc.txt2 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt
        sed -i -e "s/\s\+/\t/g" ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt
	rm -r ${i}_${REFERENCE}_fastqc
        rm ${i}_${REFERENCE}_fastqc.txt
        mv ${i}_${REFERENCE}_fastqc.html ${i}_${REFERENCE}_QC/
        mv ${i}_${REFERENCE}_fastqc.zip ${i}_${REFERENCE}_QC/
done
wait
## Contamination calculations
for i in `cat ${LIST}.txt`; do
        ONGOING=`ps aux | grep GetPileupSummaries | wc -l`
        while [ ${ONGOING} -gt ${CORES2} ]; do
                ONGOING=`ps aux | grep GetPileupSummaries | wc -l`
                echo -ne "\r${ONGOING} contamination calculation processes going on ATM, waiting for some to finish"
                sleep 1
        done
        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     ${i}: Doing contamination calculations" | tee -a ${LIST}.log
        java -Xmx10g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar GetPileupSummaries -I ${i}_${REFERENCE}.bam -V /data/Resources/References/${REFERENCE}/gnomad_biallelic.hg38.bwa.vcf.gz -O ${i}_${REFERENCE}.table.temp > /dev/null 2>&1 &
done
wait
echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Finished contamination calculations, sorting results" | tee -a ${LIST}.log
for i in `cat ${LIST}.txt`; do
        java -Xmx10g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar CalculateContamination -I ${i}_${REFERENCE}.table.temp -O ${i}_${REFERENCE}.contamination.txt > /dev/null 2>&1
        rm ${i}_${REFERENCE}.table.temp
        sed -i -e "s/whole_bam/${i}/" ${i}_${REFERENCE}.contamination.txt
        sed -i -e "s/level/sample/" ${i}_${REFERENCE}.contamination.txt
        mv ${i}_${REFERENCE}.contamination.txt ${i}_${REFERENCE}_QC/
done
wait
## Picard calculations
for i in `cat ${LIST}.txt`; do
        ONGOING=`ps aux | grep picard.jar | wc -l`
	while [ ${ONGOING} -gt ${SIMSAM} ]; do
	        ONGOING=`ps aux | grep picard.jar | wc -l`
	        echo -ne "\r${ONGOING} picard calculation processes going on ATM, waiting for some to finish"
	        sleep 1
	done
        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     ${i}: Doing picard calculations in the bed file ${BED}.bed" | tee -a ${LIST}.log
        if [ -e ${BED}.bed ] && [ ${DATATYPE} != 'WGS' ] && [ ${DATATYPE} != 'sWGS' ]; then
                java -jar /data/Resources/Software/Javas/picard.jar CollectHsMetrics R=/data/Resources/References/${REFERENCE}/${REFERENCE}.fa I=${i}_${REFERENCE}.bam O=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt BAIT_INTERVALS=${BED}.interval_list TARGET_INTERVALS=${BED}.interval_list MINIMUM_MAPPING_QUALITY=${VARMAPQ} MINIMUM_BASE_QUALITY=20 PER_TARGET_COVERAGE=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.pertarget.txt > /dev/null 2>&1 &
        else
                echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     ${i}: WARNING either bed file ${BED}.bed does not exist or the data is WGS/sWGS, so CollectHsMetrics analysis not peformed" | tee -a ${LIST}.log
        fi
        java -jar /data/Resources/Software/Javas/picard.jar CollectAlignmentSummaryMetrics R=/data/Resources/References/${REFERENCE}/${REFERENCE}.fa I=${i}_${REFERENCE}.bam O=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt > /dev/null 2>&1 &
        java -jar /data/Resources/Software/Javas/picard.jar CollectBaseDistributionByCycle CHART=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.BaseDistributionByCycle.pdf I=${i}_${REFERENCE}.bam O=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.BaseDistributionByCycle.txt > /dev/null 2>&1 &
        java -jar /data/Resources/Software/Javas/picard.jar CollectGcBiasMetrics CHART=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.GcBiasMetrics.pdf I=${i}_${REFERENCE}.bam O=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.GcBiasMetrics.txt S=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.GcBiasMetricsSummary.txt R=/data/Resources/References/${REFERENCE}/${REFERENCE}.fa > /dev/null 2>&1 &
        java -jar /data/Resources/Software/Javas/picard.jar CollectInsertSizeMetrics I=${i}_${REFERENCE}.bam O=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.InsertSizeMetrics.txt H=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.InsertSizeMetrics.pdf M=0.5 > /dev/null 2>&1 & ## have to compare the results of this with my own script for fragment size analysis
        java -jar /data/Resources/Software/Javas/picard.jar MeanQualityByCycle I=${i}_${REFERENCE}.bam O=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.MeanQualityByCycle.txt CHART=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.MeanQualityByCycle.pdf > /dev/null 2>&1 &
        java -jar /data/Resources/Software/Javas/picard.jar QualityScoreDistribution I=${i}_${REFERENCE}.bam O=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QualityScoreDistribution.txt CHART=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QualityScoreDistribution.PDF > /dev/null 2>&1 &
        if [ ${DATATYPE} = 'WGS' ] || [ ${DATATYPE} = 'sWGS' ]; then
                java -jar /data/Resources/Software/Javas/picard.jar CollectRawWgsMetrics I=${i}_${REFERENCE}.bam O=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.RawWgsMetrics.txt R=/data/Resources/References/${REFERENCE}/${REFERENCE}.fa INCLUDE_BQ_HISTOGRAM=true > /dev/null 2>&1 &
                java -jar /data/Resources/Software/Javas/picard.jar CollectWgsMetrics I=${i}_${REFERENCE}.bam O=${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt R=/data/Resources/References/${REFERENCE}/${REFERENCE}.fa > /dev/null 2>&1 &
        fi
        ## CollectTargetedPcrMetrics for amplicon panels?
        # java -jar /data/Resources/Software/Javas/picard.jar CollectTargetedPcrMetrics I=${i}_${REFERENCE}.bam O=${i}_${REFERENCE}.pcr_metrics.txt R=/data/Resources/References/${REFERENCE}/${REFERENCE}.fa AMPLICON_INTERVALS=${BED}.interval_list TARGET_INTERVALS=${BED}.interval_list 
done
wait
echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Wrapping up" | tee -a ${LIST}.log
## Annotate the pertarget file with gene and exon number
if [ -e ${BED}.bed ] && [ ${DATATYPE} != 'WGS' ] && [ ${DATATYPE} != 'sWGS' ]; then
        for i in `cat ${LIST}.txt`; do
                head -n1 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.pertarget.txt > ${i}_${REFERENCE}_QC/.head
                tail -n +2 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.pertarget.txt | awk -F '\t' '{ if ( $3 - $2 != 0 ) print }' > ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.pertarget.bed
                sort-bed ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.pertarget.bed | bedmap --echo --echo-map-id --delim '\t' - /data/Resources/BEDs/All_Exons_hg38.bed | sort -k1,1V -k2,2g > ${i}_${REFERENCE}_QC/.body
                cat ${i}_${REFERENCE}_QC/.head ${i}_${REFERENCE}_QC/.body > ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.pertarget.txt
                sed -i -e "1s/\(.\+\)/\1\tAnnotation/" ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.pertarget.txt
                rm ${i}_${REFERENCE}_QC/.head
                rm ${i}_${REFERENCE}_QC/.body
                rm ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.pertarget.bed
        done
fi
wait
## put the metrics that matter for all samples into a single file
for i in `cat ${LIST}.txt`; do
        echo -e "SAMPLE\t${i}" > ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        if [ -e ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt ]; then
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f2 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f4 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f6 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f7 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f9 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        fi
        if [ -e ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt ]; then
        head -n10 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt | tail -n4 | cut -f19 | tr "\n" "\t" | sed -e "s/\t$/\n/" | awk -F '\t' '{ print $1 "\t" $4 }' >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n10 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt | tail -n4 | cut -f20 | tr "\n" "\t" | sed -e "s/\t$/\n/" | awk -F '\t' '{ print $1 "\t" $4 }' >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n10 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt | tail -n4 | cut -f21 | tr "\n" "\t" | sed -e "s/\t$/\n/" | awk -F '\t' '{ print $1 "\t" $4 }' >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n10 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt | tail -n4 | cut -f13 | tr "\n" "\t" | sed -e "s/\t$/\n/" | awk -F '\t' '{ print $1 "\t" $4 }' >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n10 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt | tail -n4 | cut -f15 | tr "\n" "\t" | sed -e "s/\t$/\n/" | awk -F '\t' '{ print $1 "\t" $4 }' >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n10 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt | tail -n4 | cut -f22 | tr "\n" "\t" | sed -e "s/\t$/\n/" | awk -F '\t' '{ print $1 "\t" $4 }' >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n10 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt | tail -n4 | cut -f23 | tr "\n" "\t" | sed -e "s/\t$/\n/" | awk -F '\t' '{ print $1 "\t" $4 }' >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n10 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt | tail -n4 | cut -f24 | tr "\n" "\t" | sed -e "s/\t$/\n/" | awk -F '\t' '{ print $1 "\t" $4 }' >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        fi
        if [ -e ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt ]; then
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f8 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        fi
        if [ -e ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.PCRDups.txt ]; then
                if [ ${PCRDUPALG} = 'picard' ]; then
                head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.PCRDups.txt | tail -n2 | cut -f9 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
                head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.PCRDups.txt | tail -n2 | cut -f7 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
                head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.PCRDups.txt | tail -n2 | cut -f8 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
                elif [ ${PCRDUPALG} = 'samtools' ]; then
                cut -f4 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.PCRDups.txt | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
                cut -f2 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.PCRDups.txt | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
                echo -e "READ_PAIR_OPTICAL_DUPLICATES\t0" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
                fi
        fi
        if [ -e ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt ]; then
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f11 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f12 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f13 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f31 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f32 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f33 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f20 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f34 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f23 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f25 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f28 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f43 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f29 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f35 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.HsMetrics.txt | tail -n2 | cut -f53 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        fi
        if [ -e ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.contamination.txt ]; then
        cat ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.contamination.txt | cut -f2 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        cat ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.contamination.txt | cut -f3 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        fi
        if [ -e ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt ]; then
        head -n10 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt | tail -n4 | cut -f16 | tr "\n" "\t" | sed -e "s/\t$/\n/" | awk -F '\t' '{ print $1 "\t" $4 }' >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        fi
        if [ -e ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.InsertSizeMetrics.txt ]; then
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.InsertSizeMetrics.txt | tail -n2 | cut -f1 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.InsertSizeMetrics.txt | tail -n2 | cut -f3 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.InsertSizeMetrics.txt | tail -n2 | cut -f2 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.InsertSizeMetrics.txt | tail -n2 | cut -f6 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.InsertSizeMetrics.txt | tail -n2 | cut -f7 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.InsertSizeMetrics.txt | tail -n2 | cut -f9 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        fi
        if [ -e ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt ]; then
        cat ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt | cut -f2 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        cat ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt | cut -f3 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        cat ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt | cut -f4 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        cat ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt | cut -f6 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        cat ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt | cut -f7 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        cat ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt | cut -f8 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        cat ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt | cut -f8 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        cat ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt | cut -f10 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        cat ${i}_${REFERENCE}_QC/${i}_${REFERENCE}_fastqc.txt | cut -f11 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        fi
        
        if [ -e ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt ]; then
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f2 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f3 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f6 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f7 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f8 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f9 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f10 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f11 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f12 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f26 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f15 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f13 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n8 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.WgsMetrics.txt | tail -n2 | cut -f27 | tr "\n" "\t" | sed -e "s/\t$/\n/" >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        head -n10 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.AlignmentSummaryMetrics.txt | tail -n4 | cut -f2 | tr "\n" "\t" | sed -e "s/\t$/\n/" | awk -F '\t' '{ print $1 "\t" $4 }' >> ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
        fi

        awk '
        {
            for (i=1; i<=NF; i++)  {
                a[NR,i] = $i
            }
        }
        NF>p { p = NF }
        END {
            for(j=1; j<=p; j++) {
                str=a[1,j]
                for(i=2; i<=NR; i++){
                    str=str" "a[i,j];
                }
                print str
            }
        }' ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt | sed -e "s/\s\+/\t/g" > ${i}_${REFERENCE}_QC/.temp && mv ${i}_${REFERENCE}_QC/.temp ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt
done
wait
for i in `cat ${LIST}.txt`; do
        if [ -e ${LIST}_${REFERENCE}.QC.txt ]; then
                tail -n1 ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt >> ${LIST}_${REFERENCE}.QC.txt
        else
                cp ${i}_${REFERENCE}_QC/${i}_${REFERENCE}.QC.txt ${LIST}_${REFERENCE}.QC.txt
        fi
done
wait
echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Doing fragmentation analyses" | tee -a ${LIST}.log
mkdir .Temp_Frag
for SAMPLE in `cat ${LIST}.txt`; do
        ONGOING=`ps aux | grep bamToBed | wc -l`
        while [ ${ONGOING} -gt ${CORES2} ]; do
                ONGOING=`ps aux | grep bamToBed | wc -l`
                echo -ne "\r${ONGOING} frag calculation processes going on ATM, waiting for some to finish"
                sleep 1
        done
        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     ${SAMPLE}: Fragmentation analysis" | tee -a ${LIST}.log
        samtools sort -m 5G -n ${SAMPLE}_${REFERENCE}.bam 2>/dev/null | bamToBed -bedpe -i - 2>/dev/null | awk -F '\t' '{ if ( $8 > 0 ) print }' | awk -F '\t' '{ if ( $1 == $4 ) print }' > .Temp_Frag/${SAMPLE}_${REFERENCE}.bed &
done
wait
echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Wrapping up" | tee -a ${LIST}.log
if [ ! -e ${HOME}/.temp/ ]; then
        mkdir ${HOME}/.temp
fi
export TMPDIR=${HOME}/.temp
for SAMPLE in `cat ${LIST}.txt`; do
        awk -F '\t' '{ print $0 "\t" $6 - $2 }' .Temp_Frag/${SAMPLE}_${REFERENCE}.bed | sort -k1,1 -k2,2 > .Temp_Frag/${SAMPLE}_${REFERENCE}.bed.temp && mv .Temp_Frag/${SAMPLE}_${REFERENCE}.bed.temp .Temp_Frag/${SAMPLE}_${REFERENCE}.bed &
done
wait
for SAMPLE in `cat ${LIST}.txt`; do
        echo -e "Chr\tBP\tSize" > .Temp_Frag/${SAMPLE}.txt
        awk -F '\t' '{ print $1 "\t" $2 "\t" $11 }' .Temp_Frag/${SAMPLE}_${REFERENCE}.bed >> .Temp_Frag/${SAMPLE}.txt &
done
wait
echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Plotting" | tee -a ${LIST}.log
sed -i -e "s/?/0/g" ${LIST}_${REFERENCE}.QC.txt
if [ ${DATATYPE} = 'WES' ]; then
        Rscript /data/Resources/Software/ceberSUITE/IMP/R/QC_WES.r ${LIST}_${REFERENCE} 2> /dev/null
elif [ ${DATATYPE} = 'WGS' ]; then
        Rscript /data/Resources/Software/ceberSUITE/IMP/R/QC_WGS.r ${LIST}_${REFERENCE} 2> /dev/null
elif [ ${DATATYPE} = 'sWGS' ]; then
        Rscript /data/Resources/Software/ceberSUITE/IMP/R/QC_sWGS.r ${LIST}_${REFERENCE} 2> /dev/null
elif [ ${DATATYPE} = 'Capture' ]; then
        Rscript /data/Resources/Software/ceberSUITE/IMP/R/QC_Capture.r ${LIST}_${REFERENCE} 2> /dev/null
elif [ ${DATATYPE} = 'Amplicon' ]; then
        Rscript /data/Resources/Software/ceberSUITE/IMP/R/QC_Amplicon.r ${LIST}_${REFERENCE} 2> /dev/null
else
        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     R script QC: Mode not supported" | tee -a ${LIST}.log
fi
for SAMPLE in `cat ${LIST}.txt`; do
        mv ${SAMPLE}.Frag.txt ${SAMPLE}_${REFERENCE}_QC/
        mv ${SAMPLE}.Frag.png ${SAMPLE}_${REFERENCE}_QC/
done
echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Calculating fragments between 140 and 190 bp" | tee -a ${LIST}.log
for SAMPLE in `cat ${LIST}.txt`; do
        awk -F '\t' '{ if ( $11 > 140 && $11 < 190 ) print }' .Temp_Frag/${SAMPLE}_${REFERENCE}.bed | awk -F '\t' '{ print $7 }' > ${SAMPLE}_${REFERENCE}_QC/${SAMPLE}_${REFERENCE}.fragFilter.txt &
done
wait
rm -r .Temp_Frag
