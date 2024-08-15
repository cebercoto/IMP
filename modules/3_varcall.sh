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
if [ -e .AlignmentTemp_${LIST} ]; then
	rm -r .AlignmentTemp_${LIST}/
fi
mkdir .AlignmentTemp_${LIST}
CORES2=`echo "${CORES} * ${SIMSAM}" | bc -l`
## Germline Variant Calling
if [ ${SOMATIC} = 'false' ]; then
        IFS=$'\n'
        echo -e `date +%y%m%d\ %H%M%S\ %s` "     Calling variants for all individuals using the ${VARCALLALG} algorithm" | tee -a ${LIST}.log
        if [ -e ${BAMPATH}.bamlist ]; then
            for i in `cat ${LIST}.txt`; do
                    THISBAMPATH=`grep -w ${i} ${BAMPATH}.bampath | awk '{ print $2 }'`
                    if [ -e ${THISBAMPATH}/${i}_${REFERENCE}.bam ]; then
                            echo -e "${THISBAMPATH}/${i}_${REFERENCE}" >> .AlignmentTemp_${LIST}/SamplesToCall.list
                    else
                            echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Aligned BAM does not exist, connot go on" | tee -a ${LIST}.log
                    exit
                    fi
            done
        else
            for i in `cat ${LIST}.txt`; do
                if [ -e ${i}_${REFERENCE}.bam ]; then
                    echo ${i}_${REFERENCE} >> .AlignmentTemp_${LIST}/SamplesToCall.list
                else
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Aligned BAM does not exist, connot go on" | tee -a ${LIST}.log
                    exit
                fi
            done
        fi
        ## For variant calling using Unified Genotyper
        if [ ${DATATYPE} != 'sWGS' ]; then
        if [ ${VARCALLALG} = 'UG' ]; then
            if [ ${BED} = 'NULL' ]; then
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Since no bed file was provided for the variant calling with UG, using the whole genome" | tee -a ${LIST}.log
                BED=/data/Resources/References/${REFERENCE}/${REFERENCE}
            fi
            vim -c "%s/\(\S\+\)/\1\.bam/|wq" .AlignmentTemp_${LIST}/SamplesToCall.list > /dev/null 2>&1
            java -Xmx200g -jar /data/Resources/Software/Javas/GenomeAnalysisTK.jar -glm BOTH -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -T UnifiedGenotyper -D /data/Resources/References/${REFERENCE}/All.vcf.gz -o ${LIST}_${REFERENCE}.vcf -stand_call_conf 30.0 -A Coverage -A AlleleBalance --max_alternate_alleles 46 -nt ${CORES2} -I .AlignmentTemp_${LIST}/SamplesToCall.list --allow_potentially_misencoded_quality_scores -L ${BED}.bed -ip ${PAD} -dcov 1000000 -rf MappingQuality --min_mapping_quality_score ${VARMAPQ} > /dev/null 2>&1
                rm ${LIST}_${REFERENCE}.vcf.idx
                vcftools --vcf ${LIST}_${REFERENCE}.vcf --minQ 30 --min-meanDP 10 --minGQ 30 --recode --out ${LIST}_${REFERENCE}.Q.DP > /dev/null 2>&1
                vim -c "%s/\(\n\)\(\S\+\)\(\s\+\)\(\S\+\)\(\s\+\)\.\(\s\+\)/\1\2\3\4\5chr\2bp\4\6/|wq" ${LIST}_${REFERENCE}.Q.DP.recode.vcf > /dev/null 2>&1
                mv ${LIST}_${REFERENCE}.Q.DP.recode.vcf ${LIST}_${REFERENCE}.QC.vcf
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Finished calling variants" | tee -a ${LIST}.log
        fi
        ## For variant calling using Haplotype Caller, do a separate gvcf per sample, then joint variant calling
        if [ ${VARCALLALG} = 'HC' ]; then
                if [ ! -e .gvcf ]; then
                        mkdir .gvcf
                fi
                touch .AlignmentTemp_${LIST}/SamplesToCallHC.list
                for i in `cat ${LIST}.txt`; do
                        if [ -e .gvcf/${i}_${REFERENCE}.g.vcf ]; then
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Sample ${i} has already been genotyped, skipping it" | tee -a ${LIST}.log
                                echo ${i} >> .AlignmentTemp_${LIST}/SamplesToCallHC.list
                        fi
                done
                for i in `cat ${LIST}.txt | grep -v -w -f .AlignmentTemp_${LIST}/SamplesToCallHC.list`; do
                        ONGOING=`ps aux | grep HaplotypeCaller | wc -l`
                        while [ ${ONGOING} -gt ${CORES2} ]; do
                                ONGOING=`ps aux | grep HaplotypeCaller | wc -l`
                                let ONGOING2=${ONGOING}-1
                                echo -ne "\r${ONGOING2} aligments going on ATM, waiting for some to finish"
                                sleep 1
                        done
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Generating gVCF" | tee -a ${LIST}.log
                                java -Xmx20g -jar /data/Resources/Software/Javas/GenomeAnalysisTK.jar -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -T HaplotypeCaller -D /data/Resources/References/${REFERENCE}/All.vcf.gz -o .gvcf/${i}_${REFERENCE}.g.vcf -stand_call_conf 30.0 --max_alternate_alleles 46 -nct ${CORES} -I ${i}_${REFERENCE}.bam -ERC GVCF -L ${BED}.bed -ip ${PAD} -variant_index_type LINEAR -variant_index_parameter 128000 -A Coverage -A AlleleBalance --min_mapping_quality_score ${VARMAPQ} > /dev/null 2>&1 &
                                echo ${i} >> .AlignmentTemp_${LIST}/SamplesToCallHC.list
                done
                wait
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Finished generating all gVCFs" | tee -a ${LIST}.log
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Genotyping from all gVCFs" | tee -a ${LIST}.log
                sed -i -e "s/\(\S\+\)/\.gvcf\/\1_${REFERENCE}\.g\.vcf/" .AlignmentTemp_${LIST}/SamplesToCallHC.list
                java -Xmx200g -jar /data/Resources/Software/Javas/GenomeAnalysisTK.jar -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -T GenotypeGVCFs -D /data/Resources/References/${REFERENCE}/All.vcf.gz -o ${LIST}_${REFERENCE}.vcf -stand_call_conf 30.0 --max_alternate_alleles 46 -V .AlignmentTemp_${LIST}/SamplesToCallHC.list -dcov 1000000 -A Coverage -A AlleleBalance > /dev/null 2>&1
                rm ${LIST}_${REFERENCE}.vcf.idx
                rtg vcfsubset -i ${LIST}_${REFERENCE}.vcf -o ${LIST}_${REFERENCE}.vcf --remove-format PGT,PID
                gunzip -f ${LIST}_${REFERENCE}.vcf.gz
                rm ${LIST}_${REFERENCE}.vcf.gz.tbi
                vcftools --vcf ${LIST}_${REFERENCE}.vcf --minQ 30 --min-meanDP 10 --minGQ 30 --recode --out ${LIST}_${REFERENCE}.Q.DP > /dev/null 2>&1
                vim -c "%s/\(\n\)\(\S\+\)\(\s\+\)\(\S\+\)\(\s\+\)\.\(\s\+\)/\1\2\3\4\5chr\2bp\4\6/|wq" ${LIST}_${REFERENCE}.Q.DP.recode.vcf > /dev/null 2>&1
                mv ${LIST}_${REFERENCE}.Q.DP.recode.vcf ${LIST}_${REFERENCE}.QC.vcf
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Finished germline HC variant calling" | tee -a ${LIST}.log
        fi
        ## Separate each sample into one VCF was a separate module before, now incorporated here
        for i in `cat ${LIST}.txt`; do
            ONGOING=`ps aux | grep vcftools | grep ${LIST} | wc -l`
            while [ ${ONGOING} -gt ${CORES2} ]; do
                ONGOING=`ps aux | grep vcftools | grep ${LIST} | wc -l`
                let ONGOING2=${ONGOING}-1
                        echo -ne "\r${ONGOING2} aligments going on ATM, waiting for some to finish"
                sleep 1
            done
            echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Extracting indivual ${i} from the merged vcf file" | tee -a ${LIST}.log
            vcftools --vcf ${LIST}_${REFERENCE}.QC.vcf --indv ${i} --recode --out ${i}_${REFERENCE} > /dev/null 2>&1 &
        done
        wait
        echo -e `date +%y%m%d\ %H%M%S\ %s` "     Finished extracting individuals" | tee -a ${LIST}.log
        for i in `cat ${LIST}.txt`; do
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Removing missing genotypes and low GQ/DP genotypes for individual ${i}" | tee -a ${LIST}.log
                vcftools --vcf ${i}_${REFERENCE}.recode.vcf --minGQ 30 --recode --out ${i}_${REFERENCE}.GQ > /dev/null 2>&1
                vcftools --vcf ${i}_${REFERENCE}.GQ.recode.vcf --minDP 30 --recode --out ${i}_${REFERENCE}.GQ.DP > /dev/null 2>&1
                vcftools --vcf ${i}_${REFERENCE}.GQ.DP.recode.vcf --max-missing 1 --recode --out ${i}_${REFERENCE}.GQ.DP.NoMissing > /dev/null 2>&1
                mv ${i}_${REFERENCE}.recode.vcf ${i}_${REFERENCE}.vcf
                mv ${i}_${REFERENCE}.GQ.DP.NoMissing.recode.vcf ${i}_${REFERENCE}.QC.vcf
                rm ${i}_${REFERENCE}.GQ.recode.vcf
                rm ${i}_${REFERENCE}.GQ.DP.recode.vcf
        done
    fi
    ## SV calling
    if [ ${DATATYPE} != 'sWGS' ]; then
        for i in `cat ${LIST}.txt`; do
                ## Germline SV calling
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Calling Germline SVs" | tee -a ${LIST}.log
                /data/Resources/Software/delly/src/delly call -g /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -o ${i}_${REFERENCE}.delly.bcf -x /data/Resources/Software/delly/excludeTemplates/human.hg38.excl.tsv ${i}_${REFERENCE}.bam  > /dev/null 2>&1
                /data/Resources/Software/delly/src/delly merge -o sites.bcf ${i}_${REFERENCE}.delly.bcf > /dev/null 2>&1
                /data/Resources/Software/delly/src/delly call -g /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -v sites.bcf -o ${i}_${REFERENCE}.delly.bcf -x /data/Resources/Software/delly/excludeTemplates/human.hg38.excl.tsv ${i}_${REFERENCE}.bam > /dev/null 2>&1
                bcftools view ${i}_${REFERENCE}.delly.bcf | grep -v LowQual > ${i}_${REFERENCE}.delly.vcf 2>&1
                rm sites.bcf*
                rm ${i}_${REFERENCE}.delly.bcf*
                mkdir TEMP
                grep -v \# ${i}_${REFERENCE}.delly.vcf | awk -F '\t' '{ print $1 "\t" $2 "\t" $8 "\t" $10}' | sed -e "s/\(\S\+\)\s\+\S\+SVTYPE=\(\w\+\)\S\+CHR2=\(\w\+\)\S\+END=\(\w\+\)\S\+/\1\t\3\t\4\t\2/" | sed -e "s/:/\t/g" | awk -F '\t' '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $14 + $16 "\t" $15 + $17 }' | awk -F '\t' '{ print $0 "\t" $7 / ($6 + $7) "\t" NR }' > TEMP/bed1
                awk -F '\t' '{ print $1 "\t" $2 - 100 "\t" $2 + 100 "\t" $9 }' TEMP/bed1 | sort-bed - | bedmap --ec --echo --delim '\t' --echo-map-id --range 100 - /data/Resources/BEDs/All_Genes_hg38.unique.sorted2.broad.bed | awk -F '\t' '{ print $1 "_" $2 + 100 "\t" $5 "\t" $4 }' | sed -e "s/\t\t/\t\.\t/g" | sort -t$'\t' -k3,3g > TEMP/annot1
                awk -F '\t' '{ print $3 "\t" $4 - 100 "\t" $4 + 100 "\t" $9 }' TEMP/bed1 | sort-bed - | bedmap --ec --echo --delim '\t' --echo-map-id --range 100 - /data/Resources/BEDs/All_Genes_hg38.unique.sorted2.broad.bed | awk -F '\t' '{ print $1 "_" $2 + 100 "\t" $5 "\t" $4 }' | sed -e "s/\t\t/\t\.\t/g" | sort -t$'\t' -k3,3g > TEMP/annot2
                join -1 3 -2 3 TEMP/annot1 TEMP/annot2 | sed -e "s/\s\+/\t/g" | join -1 9 -2 1 TEMP/bed1 - | sed -e "s/\s\+/\t/g" | awk -F '\t' '{ print $10 "\t" $11 "\t" $12 "\t" $13 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' | sed -e "s/_/:/g" | sed -e "s/BND/TRA/g" > ${i}_${REFERENCE}.SV.bed
                if [ `cat ${i}_${REFERENCE}.SV.bed | wc -l` -gt 0 ]; then
                        sed -i -e "1i BreakPoin1\tGene1\tBreakPoint2\tGene2\tType\tReadsRef\tReadsVar\tVAF" ${i}_${REFERENCE}.SV.bed
                else
                        echo -e "BreakPoin1\tGene1\tBreakPoint2\tGene2\tType\tReadsRef\tReadsVar\tVAF" >> ${i}_${REFERENCE}.SV.bed
                fi
                rm -r TEMP
                mv ${i}_${REFERENCE}.delly.vcf ${i}_${REFERENCE}.SV.vcf
                for j in `cut -f1 ${i}_${REFERENCE}.SV.bed | sed -e "s/:/\t/" | tail -n+2`; do CHR=`echo ${j} | awk '{ print $1 }'`; BP=`echo ${j} | awk '{ print $2 }'`; TOTALREADS=`samtools view ${i}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | wc -l`; XAPAREADS=`samtools view ${i}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | grep XA: | grep -v pa: | grep "[;:]chr[0-9]*," | wc -l`; echo -e "${CHR}:${BP}\t${TOTALREADS}\t${XAPAREADS}"; done > XABP1
                for j in `cut -f3 ${i}_${REFERENCE}.SV.bed | sed -e "s/:/\t/" | tail -n+2`; do CHR=`echo ${j} | awk '{ print $1 }'`; BP=`echo ${j} | awk '{ print $2 }'`; TOTALREADS=`samtools view ${i}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | wc -l`; XAPAREADS=`samtools view ${i}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | grep XA: | grep -v pa: | grep "[;:]chr[0-9]*," | wc -l`; echo -e "${CHR}:${BP}\t${TOTALREADS}\t${XAPAREADS}"; done > XABP2
                NVARS=`cat XABP1 | wc -l`
                echo -e "Total\tXA\tXAPERC" > XA
                for j in $( seq 1 $NVARS ); do XABP1N=`head -n ${j} XABP1 | tail -n1 | awk -F '\t' '{ print $2 }'`; XABP1XA=`head -n ${j} XABP1 | tail -n1 | awk -F '\t' '{ print $3 }'`; XABP2N=`head -n ${j} XABP2 | tail -n1 | awk -F '\t' '{ print $2 }'`; XABP2XA=`head -n ${j} XABP2 | tail -n1 | awk -F '\t' '{ print $3 }'`; let TOTAL=XABP1N+XABP2N; let TOTALXA=XABP1XA+XABP2XA; XAPERC=`echo -e "${TOTAL}\t${TOTALXA}" | awk '{ print $2 / $1 }'`; echo -e "${TOTAL}\t${TOTALXA}\t${XAPERC}"; done >> XA
                paste ${i}_${REFERENCE}.SV.bed XA > XAtemp && mv XAtemp ${i}_${REFERENCE}.SV.bed
                rm XA*
        done
    fi
    ## CNV calling
    if [ ${DATATYPE} != 'sWGS' ]; then
        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: calling CNVs" | tee -a ${LIST}.log
        if [ ${DATATYPE} = 'WGS' ] || [ ${DATATYPE} = 'sWGS' ]; then
                BINLENGTH=1000
        else
                BINLENGTH=0
        fi
        mkdir .cnv_temp
        ## preprocess intervals
        echo -e `date +%y%m%d\ %H%M%S\ %s` "     Preprocessing intervals" | tee -a ${LIST}.log
        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar PreprocessIntervals -L ${BED}.interval_list -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa --bin-length ${BINLENGTH} --interval-merging-rule OVERLAPPING_ONLY -O .cnv_temp/bed.interval_list > /dev/null 2>&1 
        ######################################################## NORMAL PON for READS
        ## Count Reads on Normals
        BEDNOPATH=`echo $BED | sed -e "s/\/\S\+\///"`
        if [ -e /data/Resources/PON/${BEDNOPATH}/PON.hdf5 ]; then
                cp /data/Resources/PON/${BEDNOPATH}/PON.hdf5 .cnv_temp/${LIST}.PON.hdf5
        else
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     No option available for the PON" | tee -a ${LIST}.log
        fi
        ####################################################### END OF NORMAL PON for READS
        if [ -e .cnv_temp/${LIST}.PON.hdf5 ]; then
                ## Count Reads on Tumours
                for i in `cat ${LIST}.txt`; do
                        if [ -e ${i}_${REFERENCE}.bam ]; then
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Collecting read counts of normal sample" | tee -a ${LIST}.log
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CollectReadCounts -I ${i}_${REFERENCE}.bam -L .cnv_temp/bed.interval_list --interval-merging-rule OVERLAPPING_ONLY -O .cnv_temp/${i}_${REFERENCE}.hdf5 > /dev/null 2>&1 
                        fi
                done
                wait
                ## Denoise sample data with PON and plot
                for i in `cat ${LIST}.txt`; do
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Denoising sample using PON" | tee -a ${LIST}.log
                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar DenoiseReadCounts -I .cnv_temp/${i}_${REFERENCE}.hdf5 --count-panel-of-normals .cnv_temp/${LIST}.PON.hdf5 --standardized-copy-ratios .cnv_temp/${i}_${REFERENCE}.standardizedCR.tsv --denoised-copy-ratios .cnv_temp/${i}_${REFERENCE}.denoisedCR.tsv > /dev/null 2>&1
                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar PlotDenoisedCopyRatios --standardized-copy-ratios .cnv_temp/${i}_${REFERENCE}.standardizedCR.tsv --denoised-copy-ratios .cnv_temp/${i}_${REFERENCE}.denoisedCR.tsv --sequence-dictionary /data/Resources/References/${REFERENCE}/${REFERENCE}.dict --minimum-contig-length 46709983 --output .cnv_temp/plots --output-prefix ${i}_${REFERENCE} > /dev/null 2>&1
                done
                wait
                ## Count alleles to use in combination with the cov to call CNVs in both somatic and germline
                for i in `cat ${LIST}.txt`; do
                        ONGOING=`ps aux | grep CollectAllelicCounts | wc -l`
                        while [ ${ONGOING} -gt ${SIMSAM} ]; do
                                ONGOING=`ps aux | grep CollectAllelicCounts | wc -l`
                                echo -ne "\rProcessing ${ONGOING} samples ATM, waiting for some to finish"
                                sleep 1
                        done
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Counting alleles in normal sample" | tee -a ${LIST}.log
                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CollectAllelicCounts -L .cnv_temp/bed.interval_list -I ${i}_${REFERENCE}.bam -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -O .cnv_temp/${i}_${REFERENCE}.allelicCounts.tsv > /dev/null 2>&1 &
                done
                wait
                ## Model segments using Cov and allele counts
                for i in `cat ${LIST}.txt`; do
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Modelling segments in normal sample" | tee -a ${LIST}.log
                        java -Xmx100g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar ModelSegments --denoised-copy-ratios .cnv_temp/${i}_${REFERENCE}.denoisedCR.tsv --allelic-counts .cnv_temp/${i}_${REFERENCE}.allelicCounts.tsv --output .cnv_temp/ --output-prefix ${i}_${REFERENCE} > /dev/null 2>&1
                done               
                wait
                ## Call Segments
                for i in `cat ${LIST}.txt`; do
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Calling CNVs in normal sample" | tee -a ${LIST}.log
                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CallCopyRatioSegments --input .cnv_temp/${i}_${REFERENCE}.cr.seg --output .cnv_temp/${i}_${REFERENCE}.called.seg > /dev/null 2>&1
                        grep -v @ .cnv_temp/${i}_${REFERENCE}.called.seg | awk -F '\t' '{ if ( $5 >= 0.2 || $5 <= -0.2 ) print }' > .cnv_temp/${i}_${REFERENCE}.CNV.bed
                done
                wait
                ## Plot results
                for i in `cat ${LIST}.txt`; do
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Plotting CNVs" | tee -a ${LIST}.log
                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar PlotModeledSegments --denoised-copy-ratios .cnv_temp/${i}_${REFERENCE}.denoisedCR.tsv --allelic-counts .cnv_temp/${i}_${REFERENCE}.hets.tsv -segments .cnv_temp/${i}_${REFERENCE}.modelFinal.seg --sequence-dictionary /data/Resources/References/${REFERENCE}/${REFERENCE}.dict --minimum-contig-length 46709983 --output .cnv_temp/plots --output-prefix ${i}_${REFERENCE} > /dev/null 2>&1
                        cp .cnv_temp/plots/${i}_${REFERENCE}.modeled.png ${i}_${REFERENCE}.CNV.png
                done
                wait
                for i in `cat ${LIST}.txt`; do
                        head -n1 .cnv_temp/${i}_${REFERENCE}.CNV.bed > .cnv_temp/head
                        sed -i -e "s/\(.\+\)/\1\tGENES/" .cnv_temp/head
                        tail -n+2 .cnv_temp/${i}_${REFERENCE}.CNV.bed | sort-bed - | bedmap --ec --echo --delim '\t' --echo-map-id --range 100 - /data/Resources/BEDs/All_Genes_hg38.unique.sorted2.broad.bed > .cnv_temp/tail
                        cat .cnv_temp/head .cnv_temp/tail > ${i}_${REFERENCE}.CNV.bed
                done
                wait
                rm -r .cnv_temp/      
        else
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     No option available for CNV PON. Not calling CNVs." | tee -a ${LIST}.log
        fi
    fi
    echo -e `date +%y%m%d\ %H%M%S\ %s` "     Finished calling variants" | tee -a ${LIST}.log
## The somatic Variant calling pipeline
elif [ ${SOMATIC} = 'true' ]; then
        echo -e `date +%y%m%d\ %H%M%S\ %s` "     Calling variants for all individuals using the Somatic (Mutect2/UG|HC) Pipeline" | tee -a ${LIST}.log
        ######################### have to implement here the bam list thing for full path of bams #######################
        IFS=$'\n'
        if [ ${DATATYPE} != 'sWGS' ]; then
            if [ ${VARCALLALG} != 'None' ]; then
                ## Prepare germline calls
                touch .AlignmentTemp_${LIST}/SamplesToCallGL.list
                for i in `cat ${SOMATICFILE}.txt | awk -F '\t' '{ print $1}' | sort -u`; do
                        if [ -e ${i}_${REFERENCE}.bam ]; then
                            echo ${i} >> .AlignmentTemp_${LIST}/SamplesToCallGL.list
                        fi
                done
                if [ ${BED} = 'NULL' ]; then
                    echo -e `date +%y%m%d\ %H%M%S\ %s` "     Since no bed file was provided for the variant calling with UG, using the whole genome" | tee -a ${LIST}.log
                    BED=/data/Resources/References/${REFERENCE}/${REFERENCE}
                fi
                ## UG germline calls
                if [ ${VARCALLALG} = 'UG' ]; then
                        if [ `cat .AlignmentTemp_${LIST}/SamplesToCallGL.list | wc -l` -gt 0 ]; then
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Calling variants for all individuals using the UG algorithm" | tee -a ${LIST}.log
                                cp .AlignmentTemp_${LIST}/SamplesToCallGL.list .AlignmentTemp_${LIST}/SamplesToCallUG.list
                                vim -c "%s/\(\S\+\)/\1_${REFERENCE}\.bam/|wq" .AlignmentTemp_${LIST}/SamplesToCallUG.list > /dev/null 2>&1
                                java -Xmx200g -jar /data/Resources/Software/Javas/GenomeAnalysisTK.jar -glm BOTH -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -T UnifiedGenotyper -D /data/Resources/References/${REFERENCE}/All.vcf.gz -o ${LIST}_${REFERENCE}.vcf -stand_call_conf 30.0 -A Coverage -A AlleleBalance --max_alternate_alleles 46 -nt ${CORES2} -I .AlignmentTemp_${LIST}/SamplesToCallUG.list --allow_potentially_misencoded_quality_scores -L ${BED}.bed -ip ${PAD} -dcov 1000000 -rf MappingQuality --min_mapping_quality_score ${VARMAPQ} > /dev/null 2>&1
                                rm ${LIST}_${REFERENCE}.vcf.idx
                                vcftools --vcf ${LIST}_${REFERENCE}.vcf --minQ 30 --min-meanDP 10 --minGQ 30 --recode --out ${LIST}_${REFERENCE}.Q.DP > /dev/null 2>&1
                                vim -c "%s/\(\n\)\(\S\+\)\(\s\+\)\(\S\+\)\(\s\+\)\.\(\s\+\)/\1\2\3\4\5chr\2bp\4\6/|wq" ${LIST}_${REFERENCE}.Q.DP.recode.vcf > /dev/null 2>&1
                                mv ${LIST}_${REFERENCE}.Q.DP.recode.vcf ${LIST}_${REFERENCE}.QC.vcf
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Finished germline UG variant calling" | tee -a ${LIST}.log
                        fi
                fi
                ## HC germline calls
                if [ ${VARCALLALG} = 'HC' ]; then
                        if [ `cat .AlignmentTemp_${LIST}/SamplesToCallGL.list | wc -l` -gt 0 ]; then
                                if [ ! -e .gvcf ]; then
                                        mkdir .gvcf
                                fi
                                touch .AlignmentTemp_${LIST}/SamplesToCallHC.list
                                for i in `cat .AlignmentTemp_${LIST}/SamplesToCallGL.list`; do
                                        if [ -e .gvcf/${i}_${REFERENCE}.g.vcf ]; then
                                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Sample ${i} has already been genotyped, skipping it" | tee -a ${LIST}.log
                                                echo ${i} >> .AlignmentTemp_${LIST}/SamplesToCallHC.list
                                        fi
                                done
                                for i in `cat .AlignmentTemp_${LIST}/SamplesToCallGL.list | grep -v -w -f .AlignmentTemp_${LIST}/SamplesToCallHC.list`; do
                                        ONGOING=`ps aux | grep HaplotypeCaller | wc -l`
                                        while [ ${ONGOING} -gt ${CORES2} ]; do
                                                ONGOING=`ps aux | grep HaplotypeCaller | wc -l`
                                                let ONGOING2=${ONGOING}-1
                                                echo -ne "\r${ONGOING2} variant calls going on ATM, waiting for some to finish"
                                                sleep 1
                                        done
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Generating gVCF" | tee -a ${LIST}.log
                                        java -Xmx20g -jar /data/Resources/Software/Javas/GenomeAnalysisTK.jar -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -T HaplotypeCaller -D /data/Resources/References/${REFERENCE}/All.vcf.gz -o .gvcf/${i}_${REFERENCE}.g.vcf -stand_call_conf 30.0 --max_alternate_alleles 46 -nct ${CORES} -I ${i}_${REFERENCE}.bam -ERC GVCF -L ${BED}.bed -ip ${PAD} -variant_index_type LINEAR -variant_index_parameter 128000 -A Coverage -A AlleleBalance --min_mapping_quality_score ${VARMAPQ} > /dev/null 2>&1 &
                                        echo ${i} >> .AlignmentTemp_${LIST}/SamplesToCallHC.list
                                        sleep 2
                                done
                                wait
                                sed -i -e "s/\(\S\+\)/.gvcf\/\1_${REFERENCE}.g.vcf/" .AlignmentTemp_${LIST}/SamplesToCallHC.list
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Finished generating all gVCFs" | tee -a ${LIST}.log
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Genotyping from all gVCFs" | tee -a ${LIST}.log
                                java -Xmx200g -jar /data/Resources/Software/Javas/GenomeAnalysisTK.jar -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -T GenotypeGVCFs -D /data/Resources/References/${REFERENCE}/All.vcf.gz -o ${LIST}_${REFERENCE}.vcf -stand_call_conf 30.0 --max_alternate_alleles 46 -nt ${CORES} -V .AlignmentTemp_${LIST}/SamplesToCallHC.list -dcov 1000000 -A Coverage -A AlleleBalance > /dev/null 2>&1
                                #rm ${LIST}_${REFERENCE}.vcf.idx
                                rtg vcfsubset -i ${LIST}_${REFERENCE}.vcf -o ${LIST}_${REFERENCE}.vcf --remove-format PGT,PID
                                gunzip -f ${LIST}_${REFERENCE}.vcf.gz
                                vcftools --vcf ${LIST}_${REFERENCE}.vcf --minQ 30 --min-meanDP 10 --minGQ 30 --recode --out ${LIST}_${REFERENCE}.Q.DP > /dev/null 2>&1
                                vim -c "%s/\(\n\)\(\S\+\)\(\s\+\)\(\S\+\)\(\s\+\)\.\(\s\+\)/\1\2\3\4\5chr\2bp\4\6/|wq" ${LIST}_${REFERENCE}.Q.DP.recode.vcf > /dev/null 2>&1
                                mv ${LIST}_${REFERENCE}.Q.DP.recode.vcf ${LIST}_${REFERENCE}.QC.vcf
                                rm *.tbi
                            echo -e `date +%y%m%d\ %H%M%S\ %s` "     Finished germline HC variant calling" | tee -a ${LIST}.log
                        fi
                fi
                ## Separate and filter individual VCFs
                if [ `cat .AlignmentTemp_${LIST}/SamplesToCallGL.list | wc -l` -gt 0 ]; then
                    for i in `cat .AlignmentTemp_${LIST}/SamplesToCallGL.list`; do
                        ONGOING=`ps aux | grep vcftools | grep ${LIST} | wc -l`
                        while [ ${ONGOING} -gt ${CORES2} ]; do
                                ONGOING=`ps aux | grep vcftools | grep ${LIST} | wc -l`
                                echo -ne "\r${ONGOING} Individual extractions going on ATM, waiting for some to finish"
                                sleep 1
                        done
                        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Extracting indivual ${i} from the merged vcf file" | tee -a ${LIST}.log
                        vcftools --vcf ${LIST}_${REFERENCE}.QC.vcf --indv ${i} --recode --out ${i}_${REFERENCE} > /dev/null 2>&1 &
                    done
                    wait
                    echo -e `date +%y%m%d\ %H%M%S\ %s` "     Finished extracting individuals" | tee -a ${LIST}.log
                fi
                for i in `cat .AlignmentTemp_${LIST}/SamplesToCallGL.list`; do
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     Removing missing genotypes and low GQ/DP genotypes for individual ${i}" | tee -a ${LIST}.log
                        vcftools --vcf ${i}_${REFERENCE}.recode.vcf --minGQ 30 --recode --out ${i}_${REFERENCE}.GQ > /dev/null 2>&1
                        vcftools --vcf ${i}_${REFERENCE}.GQ.recode.vcf --minDP 30 --recode --out ${i}_${REFERENCE}.GQ.DP > /dev/null 2>&1
                        vcftools --vcf ${i}_${REFERENCE}.GQ.DP.recode.vcf --max-missing 1 --recode --out ${i}_${REFERENCE}.GQ.DP.NoMissing > /dev/null 2>&1
                        mv ${i}_${REFERENCE}.recode.vcf ${i}_${REFERENCE}.vcf
                        mv ${i}_${REFERENCE}.GQ.DP.NoMissing.recode.vcf ${i}_${REFERENCE}.QC.vcf
                        rm ${i}_${REFERENCE}.GQ.recode.vcf
                        rm ${i}_${REFERENCE}.GQ.DP.recode.vcf
                done
                ## BEGINNNING FOR SOMATIC SIDE OF THINGS: Creation of PON
                if [ `cat .AlignmentTemp_${LIST}/SamplesToCallGL.list | wc -l` -gt 0 ]; then
                        if [ ! -e PON.vcf ]; then
                                for i in `cat ${SOMATICFILE}.txt`; do
                                        ONGOING=`ps aux | grep Mutect2 | wc -l`
                                        while [ ${ONGOING} -gt ${CORES2} ]; do
                                                ONGOING=`ps aux | grep Mutect2 | wc -l`
                                                echo -ne "\rProcessing ${ONGOING} samples ATM, waiting for some to finish"
                                                sleep 1
                                        done
                                        NORMAL=`echo ${i} | awk '{ print $1 }'`
                                        TUMOUR=`echo ${i} | awk '{ print $2 }'`
                                        if [ ! -e ${NORMAL}_${REFERENCE}_forPON.vcf ]; then
                                                if [ -e ${NORMAL}_${REFERENCE}.bam ]; then
                                                        echo ${NORMAL}_${REFERENCE}_forPON.vcf >> PON.list
                                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${NORMAL}: Calling variants for PON creation (tumour pair is ${TUMOUR})" | tee -a ${LIST}.log
                                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar Mutect2 -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -I ${NORMAL}_${REFERENCE}.bam -tumor ${NORMAL} -O ${NORMAL}_${REFERENCE}_forPON.vcf -L ${BED}.bed -ip ${PAD} > /dev/null 2>&1 &
                                                        touch ${NORMAL}_${REFERENCE}_forPON.vcf
                                                else
                                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${TUMOUR}: No normal provided" | tee -a ${LIST}.log
                                                fi
                                        fi
                                        sleep 4
                                done
                                wait
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Creating PON panel" | tee -a ${LIST}.log
                                if [ -e PON.list ]; then
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar CreateSomaticPanelOfNormals -vcfs PON.list -O PON.vcf > /dev/null 2>&1
                                        rm *forPON.vcf
                                        rm PON.list
                                fi
                        fi
                else
                        BEDNOPATH=`echo $BED | sed -e "s/\/\S\+\///"`
                        if [ -e /data/Resources/PON/${BEDNOPATH}/PON.vcf ]; then
                                cp /data/Resources/PON/${BEDNOPATH}/PON.vcf PON.vcf
                                cp /data/Resources/PON/${BEDNOPATH}/PON.vcf.idx PON.vcf.idx
                        else
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     No option available for the PON" | tee -a ${LIST}.log
                        fi
                fi
                ## MuTect2 Variant Calling for each pair (and if there's no germline also for the somatic without using the germline)
                for i in `cat ${SOMATICFILE}.txt`; do
                        ONGOING=`ps aux | grep Mutect2 | wc -l`
                        while [ ${ONGOING} -gt ${CORES2} ]; do
                                ONGOING=`ps aux | grep Mutect2 | wc -l`
                                echo -ne "\rProcessing ${ONGOING} samples ATM, waiting for some to finish"
                                sleep 1
                        done
                        NORMAL=`echo ${i} | awk '{ print $1 }'`
                        TUMOUR=`echo ${i} | awk '{ print $2 }'`
                        if [ -e ${NORMAL}_${REFERENCE}.bam ]; then
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${TUMOUR}: Calling variants using normal paired ${NORMAL}" | tee -a ${LIST}.log
                                if [ -e PON.vcf ]; then
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar Mutect2 -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -I ${TUMOUR}_${REFERENCE}.bam -tumor ${TUMOUR} -I ${NORMAL}_${REFERENCE}.bam -normal ${NORMAL} --germline-resource /data/Resources/References/${REFERENCE}/gnomad.hg38.bwa.vcf.gz --panel-of-normals PON.vcf -O ${TUMOUR}_${REFERENCE}_somatic.vcf -L ${BED}.bed -ip ${PAD} --max-reads-per-alignment-start 0 > /dev/null 2>&1 &
                                else
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar Mutect2 -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -I ${TUMOUR}_${REFERENCE}.bam -tumor ${TUMOUR} -I ${NORMAL}_${REFERENCE}.bam -normal ${NORMAL} --germline-resource /data/Resources/References/${REFERENCE}/gnomad.hg38.bwa.vcf.gz -O ${TUMOUR}_${REFERENCE}_somatic.vcf -L ${BED}.bed -ip ${PAD} --max-reads-per-alignment-start 0 > /dev/null 2>&1 &
                                fi
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar Mutect2 -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -I ${TUMOUR}_${REFERENCE}.bam -tumor ${TUMOUR} -O ${TUMOUR}_${REFERENCE}_somaticNOGLF.vcf -L ${BED}.bed -ip ${PAD} --max-reads-per-alignment-start 0 > /dev/null 2>&1 &
                        else
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${TUMOUR}: Calling variants without a normal" | tee -a ${LIST}.log
                                if [ -e PON.vcf ]; then
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar Mutect2 -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -I ${TUMOUR}_${REFERENCE}.bam -tumor ${TUMOUR} --germline-resource /data/Resources/References/${REFERENCE}/gnomad.hg38.bwa.vcf.gz --panel-of-normals PON.vcf -O ${TUMOUR}_${REFERENCE}_somatic.vcf -L ${BED}.bed -ip ${PAD} --max-reads-per-alignment-start 0 > /dev/null 2>&1 &
                                else
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar Mutect2 -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -I ${TUMOUR}_${REFERENCE}.bam -tumor ${TUMOUR} --germline-resource /data/Resources/References/${REFERENCE}/gnomad.hg38.bwa.vcf.gz -O ${TUMOUR}_${REFERENCE}_somatic.vcf -L ${BED}.bed -ip ${PAD} --max-reads-per-alignment-start 0 > /dev/null 2>&1 &
                                fi
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar Mutect2 -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -I ${TUMOUR}_${REFERENCE}.bam -tumor ${TUMOUR} -O ${TUMOUR}_${REFERENCE}_somaticNOGLF.vcf -L ${BED}.bed -ip ${PAD} --max-reads-per-alignment-start 0 > /dev/null 2>&1 &
                        fi
                done
                wait
                for i in `cat ${SOMATICFILE}.txt`; do
                        NORMAL=`echo ${i} | awk '{ print $1 }'`
                        TUMOUR=`echo ${i} | awk '{ print $2 }'`
                        rtg vcfsubset -i ${TUMOUR}_${REFERENCE}_somatic.vcf -o ${TUMOUR}_${REFERENCE}_somatic.vcf --remove-format PGT,PID
                        rtg vcfsubset -i ${TUMOUR}_${REFERENCE}_somaticNOGLF.vcf -o ${TUMOUR}_${REFERENCE}_somaticNOGLF.vcf --remove-format PGT,PID
                        gunzip -f ${TUMOUR}_${REFERENCE}_somatic.vcf.gz
                        gunzip -f ${TUMOUR}_${REFERENCE}_somaticNOGLF.vcf.gz
                done
                wait
                rm *.tbi
                ## Apply filters to the VCF
                for i in `cat ${LIST}.txt`; do
                        ONGOING=`ps aux | grep GetPileupSummaries | wc -l`
                        while [ ${ONGOING} -gt ${CORES2} ]; do
                                ONGOING=`ps aux | grep GetPileupSummaries | wc -l`
                                echo -ne "\rProcessing ${ONGOING} samples ATM, waiting for some to finish"
                                sleep 1
                        done
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Calculating stats for filters" | tee -a ${LIST}.log
                        java -Xmx10g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar GetPileupSummaries -I ${i}_${REFERENCE}.bam -V /data/Resources/References/${REFERENCE}/gnomad_biallelic.hg38.bwa.vcf.gz -O ${i}_${REFERENCE}.table > /dev/null 2>&1 &
                done
                wait
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Applying filters" | tee -a ${LIST}.log
                for i in `cat ${SOMATICFILE}.txt`; do
                        NORMAL=`echo ${i} | awk '{ print $1 }'`
                        TUMOUR=`echo ${i} | awk '{ print $2 }'`
                        if [ -e ${NORMAL}_${REFERENCE}.bam ]; then
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar CalculateContamination -I ${TUMOUR}_${REFERENCE}.table -O ${TUMOUR}_${REFERENCE}.contamination.table -matched ${NORMAL}_${REFERENCE}.table > /dev/null 2>&1
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar FilterMutectCalls -V ${TUMOUR}_${REFERENCE}_somatic.vcf --contamination-table ${TUMOUR}_${REFERENCE}.contamination.table -O ${TUMOUR}_${REFERENCE}_somatic.QC.vcf > /dev/null 2>&1
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar CalculateContamination -I ${TUMOUR}_${REFERENCE}.table -O ${TUMOUR}_${REFERENCE}.contamination.table -matched ${NORMAL}_${REFERENCE}.table > /dev/null 2>&1
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar FilterMutectCalls -V ${TUMOUR}_${REFERENCE}_somaticNOGLF.vcf --contamination-table ${TUMOUR}_${REFERENCE}.contamination.table -O ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.vcf > /dev/null 2>&1
                        else
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar CalculateContamination -I ${TUMOUR}_${REFERENCE}.table -O ${TUMOUR}_${REFERENCE}.contamination.table > /dev/null 2>&1
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar FilterMutectCalls -V ${TUMOUR}_${REFERENCE}_somatic.vcf --contamination-table ${TUMOUR}_${REFERENCE}.contamination.table -O ${TUMOUR}_${REFERENCE}_somatic.QC.vcf > /dev/null 2>&1
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar CalculateContamination -I ${TUMOUR}_${REFERENCE}.table -O ${TUMOUR}_${REFERENCE}.contamination.table > /dev/null 2>&1
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.0.3.0-local.jar FilterMutectCalls -V ${TUMOUR}_${REFERENCE}_somaticNOGLF.vcf --contamination-table ${TUMOUR}_${REFERENCE}.contamination.table -O ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.vcf > /dev/null 2>&1
                        fi
                done
                wait
                rm *.table
                rm *.idx
                for i in `cat ${SOMATICFILE}.txt`; do
                        NORMAL=`echo ${i} | awk '{ print $1 }'`
                        TUMOUR=`echo ${i} | awk '{ print $2 }'`
                        if [ -e ${NORMAL}_${REFERENCE}.bam ]; then
                                bcftools view -Ov -s ${TUMOUR} ${TUMOUR}_${REFERENCE}_somatic.vcf > ${TUMOUR}_${REFERENCE}_somatic.recode.vcf
                                mv ${TUMOUR}_${REFERENCE}_somatic.recode.vcf ${TUMOUR}_${REFERENCE}_somatic.vcf
                                bcftools view -Ov -s ${TUMOUR} ${TUMOUR}_${REFERENCE}_somatic.QC.vcf > ${TUMOUR}_${REFERENCE}_somatic.QC.recode.vcf
                                mv ${TUMOUR}_${REFERENCE}_somatic.QC.recode.vcf ${TUMOUR}_${REFERENCE}_somatic.QC.vcf
                                bcftools view -Ov -s ${TUMOUR} ${TUMOUR}_${REFERENCE}_somaticNOGLF.vcf > ${TUMOUR}_${REFERENCE}_somaticNOGLF.recode.vcf
                                mv ${TUMOUR}_${REFERENCE}_somaticNOGLF.recode.vcf ${TUMOUR}_${REFERENCE}_somaticNOGLF.vcf
                                bcftools view -Ov -s ${TUMOUR} ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.vcf > ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.recode.vcf
                                mv ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.recode.vcf ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.vcf
                        fi
                done
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Wrapping up" | tee -a ${LIST}.log
                for i in `cat ${SOMATICFILE}.txt`; do
                        NORMAL=`echo ${i} | awk '{ print $1 }'`
                        TUMOUR=`echo ${i} | awk '{ print $2 }'`
                        vim -c "%s/\(\n\)\(\S\+\)\(\s\+\)\(\S\+\)\(\s\+\)\.\(\s\+\)/\1\2\3\4\5chr\2bp\4\6/e|wq" ${TUMOUR}_${REFERENCE}_somatic.vcf > /dev/null 2>&1
                        vim -c "%s/\(\n\)\(\S\+\)\(\s\+\)\(\S\+\)\(\s\+\)\.\(\s\+\)/\1\2\3\4\5chr\2bp\4\6/e|wq" ${TUMOUR}_${REFERENCE}_somatic.QC.vcf > /dev/null 2>&1
                        vim -c "%s/\(\n\)\(\S\+\)\(\s\+\)\(\S\+\)\(\s\+\)\.\(\s\+\)/\1\2\3\4\5chr\2bp\4\6/e|wq" ${TUMOUR}_${REFERENCE}_somaticNOGLF.vcf > /dev/null 2>&1
                        vim -c "%s/\(\n\)\(\S\+\)\(\s\+\)\(\S\+\)\(\s\+\)\.\(\s\+\)/\1\2\3\4\5chr\2bp\4\6/e|wq" ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.vcf > /dev/null 2>&1
                done
                wait
            fi
        fi
        ## SV calling
        if [ ${DATATYPE} != 'sWGS' ]; then
            for i in `cat ${SOMATICFILE}.txt`; do
                    NORMAL=`echo ${i} | awk '{ print $1 }'`
                    TUMOUR=`echo ${i} | awk '{ print $2 }'`
                    if [ -e ${NORMAL}_${REFERENCE}.bam ]; then
                            ## Germline SV calling
                            if [ ! -e ${NORMAL}_${REFERENCE}.SV.bed ]; then
                                    echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Calling Germline SVs" | tee -a ${LIST}.log
                                    /data/Resources/Software/delly/src/delly call -g /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -o ${NORMAL}_${REFERENCE}.delly.bcf -x /data/Resources/Software/delly/excludeTemplates/human.hg38.excl.tsv ${NORMAL}_${REFERENCE}.bam  > /dev/null 2>&1
                                    /data/Resources/Software/delly/src/delly merge -o sites.bcf ${NORMAL}_${REFERENCE}.delly.bcf > /dev/null 2>&1
                                    /data/Resources/Software/delly/src/delly call -g /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -v sites.bcf -o ${NORMAL}_${REFERENCE}.delly.bcf -x /data/Resources/Software/delly/excludeTemplates/human.hg38.excl.tsv ${NORMAL}_${REFERENCE}.bam > /dev/null 2>&1
                                    bcftools view ${NORMAL}_${REFERENCE}.delly.bcf | grep -v LowQual > ${NORMAL}_${REFERENCE}.delly.vcf 2>&1
                                    rm sites.bcf*
                                    rm ${NORMAL}_${REFERENCE}.delly.bcf*
                                    mkdir TEMP
                                    grep -v \# ${NORMAL}_${REFERENCE}.delly.vcf | awk -F '\t' '{ print $1 "\t" $2 "\t" $8 "\t" $10}' | sed -e "s/\(\S\+\)\s\+\S\+SVTYPE=\(\w\+\)\S\+CHR2=\(\w\+\)\S\+END=\(\w\+\)\S\+/\1\t\3\t\4\t\2/" | sed -e "s/:/\t/g" | awk -F '\t' '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $14 + $16 "\t" $15 + $17 }' | awk -F '\t' '{ print $0 "\t" $7 / ($6 + $7) "\t" NR }' > TEMP/bed1
                                    awk -F '\t' '{ print $1 "\t" $2 - 100 "\t" $2 + 100 "\t" $9 }' TEMP/bed1 | sort-bed - | bedmap --ec --echo --delim '\t' --echo-map-id --range 100 - /data/Resources/BEDs/All_Genes_hg38.unique.sorted2.broad.bed | awk -F '\t' '{ print $1 "_" $2 + 100 "\t" $5 "\t" $4 }' | sed -e "s/\t\t/\t\.\t/g" | sort -t$'\t' -k3,3g > TEMP/annot1
                                    awk -F '\t' '{ print $3 "\t" $4 - 100 "\t" $4 + 100 "\t" $9 }' TEMP/bed1 | sort-bed - | bedmap --ec --echo --delim '\t' --echo-map-id --range 100 - /data/Resources/BEDs/All_Genes_hg38.unique.sorted2.broad.bed | awk -F '\t' '{ print $1 "_" $2 + 100 "\t" $5 "\t" $4 }' | sed -e "s/\t\t/\t\.\t/g" | sort -t$'\t' -k3,3g > TEMP/annot2
                                    join -1 3 -2 3 TEMP/annot1 TEMP/annot2 | sed -e "s/\s\+/\t/g" | join -1 9 -2 1 TEMP/bed1 - | sed -e "s/\s\+/\t/g" | awk -F '\t' '{ print $10 "\t" $11 "\t" $12 "\t" $13 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' | sed -e "s/_/:/g" | sed -e "s/BND/TRA/g" > ${NORMAL}_${REFERENCE}.SV.bed
                                    if [ `cat ${NORMAL}_${REFERENCE}.SV.bed | wc -l` -gt 0 ]; then
                                            sed -i -e "1i BreakPoin1\tGene1\tBreakPoint2\tGene2\tType\tReadsRef\tReadsVar\tVAF" ${NORMAL}_${REFERENCE}.SV.bed
                                    else
                                            echo -e "BreakPoin1\tGene1\tBreakPoint2\tGene2\tType\tReadsRef\tReadsVar\tVAF" >> ${NORMAL}_${REFERENCE}.SV.bed
                                    fi
                                    rm -r TEMP
                                    mv ${NORMAL}_${REFERENCE}.delly.vcf ${NORMAL}_${REFERENCE}.SV.vcf
                                    for j in `cut -f1 ${NORMAL}_${REFERENCE}.SV.bed | sed -e "s/:/\t/" | tail -n+2`; do CHR=`echo ${j} | awk '{ print $1 }'`; BP=`echo ${j} | awk '{ print $2 }'`; TOTALREADS=`samtools view ${NORMAL}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | wc -l`; XAPAREADS=`samtools view ${NORMAL}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | grep XA: | grep -v pa: | grep "[;:]chr[0-9]*," | wc -l`; echo -e "${CHR}:${BP}\t${TOTALREADS}\t${XAPAREADS}"; done > XABP1
                                    for j in `cut -f3 ${NORMAL}_${REFERENCE}.SV.bed | sed -e "s/:/\t/" | tail -n+2`; do CHR=`echo ${j} | awk '{ print $1 }'`; BP=`echo ${j} | awk '{ print $2 }'`; TOTALREADS=`samtools view ${NORMAL}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | wc -l`; XAPAREADS=`samtools view ${NORMAL}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | grep XA: | grep -v pa: | grep "[;:]chr[0-9]*," | wc -l`; echo -e "${CHR}:${BP}\t${TOTALREADS}\t${XAPAREADS}"; done > XABP2
                                    NVARS=`cat XABP1 | wc -l`
                                    echo -e "Total\tXA\tXAPERC" > XA
                                    for j in $( seq 1 $NVARS ); do XABP1N=`head -n ${j} XABP1 | tail -n1 | awk -F '\t' '{ print $2 }'`; XABP1XA=`head -n ${j} XABP1 | tail -n1 | awk -F '\t' '{ print $3 }'`; XABP2N=`head -n ${j} XABP2 | tail -n1 | awk -F '\t' '{ print $2 }'`; XABP2XA=`head -n ${j} XABP2 | tail -n1 | awk -F '\t' '{ print $3 }'`; let TOTAL=XABP1N+XABP2N; let TOTALXA=XABP1XA+XABP2XA; XAPERC=`echo -e "${TOTAL}\t${TOTALXA}" | awk '{ print $2 / $1 }'`; echo -e "${TOTAL}\t${TOTALXA}\t${XAPERC}"; done >> XA
                                    paste ${NORMAL}_${REFERENCE}.SV.bed XA > XAtemp && mv XAtemp ${NORMAL}_${REFERENCE}.SV.bed
                                    rm XA*
                            fi
                            ## Somatic SV calling
                            echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Calling Tumour SVs" | tee -a ${LIST}.log
                            echo -e "${TUMOUR}\ttumor" >> ${TUMOUR}.pair.txt
                            echo -e "${NORMAL}\tcontrol" >> ${TUMOUR}.pair.txt
                            /data/Resources/Software/delly/src/delly call -x /data/Resources/Software/delly/excludeTemplates/human.hg38.excl.tsv -o ${TUMOUR}_${REFERENCE}.delly.bcf -g /data/Resources/References/${REFERENCE}/${REFERENCE}.fa ${TUMOUR}_${REFERENCE}.bam ${NORMAL}_${REFERENCE}.bam > /dev/null 2>&1
                            /data/Resources/Software/delly/src/delly filter -f somatic -o ${TUMOUR}_${REFERENCE}.delly.pre.bcf -s ${TUMOUR}.pair.txt ${TUMOUR}_${REFERENCE}.delly.bcf > /dev/null 2>&1
                            /data/Resources/Software/delly/src/delly call -g /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -v ${TUMOUR}_${REFERENCE}.delly.pre.bcf -o ${TUMOUR}_${REFERENCE}.delly.bcf -x /data/Resources/Software/delly/excludeTemplates/human.hg38.excl.tsv ${TUMOUR}_${REFERENCE}.bam ${NORMAL}_${REFERENCE}.bam > /dev/null 2>&1
                            /data/Resources/Software/delly/src/delly filter -f somatic -o ${TUMOUR}_${REFERENCE}.delly.filter.bcf -s ${TUMOUR}.pair.txt ${TUMOUR}_${REFERENCE}.delly.bcf > /dev/null 2>&1
                            bcftools view ${TUMOUR}_${REFERENCE}.delly.filter.bcf | grep -v LowQual > ${TUMOUR}_${REFERENCE}.delly.vcf 2>&1
                            rm *.bcf
                            rm *.csi
                            rm ${TUMOUR}.pair.txt
                            mkdir TEMP
                            grep -v \# ${TUMOUR}_${REFERENCE}.delly.vcf | awk -F '\t' '{ print $1 "\t" $2 "\t" $8 "\t" $10}' | sed -e "s/\(\S\+\)\s\+\S\+SVTYPE=\(\w\+\)\S\+CHR2=\(\w\+\)\S\+END=\(\w\+\)\S\+/\1\t\3\t\4\t\2/" | sed -e "s/:/\t/g" | awk -F '\t' '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $14 + $16 "\t" $15 + $17 }' | awk -F '\t' '{ print $0 "\t" $7 / ($6 + $7) "\t" NR }' > TEMP/bed1
                            awk -F '\t' '{ print $1 "\t" $2 - 100 "\t" $2 + 100 "\t" $9 }' TEMP/bed1 | sort-bed - | bedmap --ec --echo --delim '\t' --echo-map-id --range 100 - /data/Resources/BEDs/All_Genes_hg38.unique.sorted2.broad.bed | awk -F '\t' '{ print $1 "_" $2 + 100 "\t" $5 "\t" $4 }' | sed -e "s/\t\t/\t\.\t/g" | sort -t$'\t' -k3,3g > TEMP/annot1
                            awk -F '\t' '{ print $3 "\t" $4 - 100 "\t" $4 + 100 "\t" $9 }' TEMP/bed1 | sort-bed - | bedmap --ec --echo --delim '\t' --echo-map-id --range 100 - /data/Resources/BEDs/All_Genes_hg38.unique.sorted2.broad.bed | awk -F '\t' '{ print $1 "_" $2 + 100 "\t" $5 "\t" $4 }' | sed -e "s/\t\t/\t\.\t/g" | sort -t$'\t' -k3,3g > TEMP/annot2
                            join -1 3 -2 3 TEMP/annot1 TEMP/annot2 | sed -e "s/\s\+/\t/g" | join -1 9 -2 1 TEMP/bed1 - | sed -e "s/\s\+/\t/g" | awk -F '\t' '{ print $10 "\t" $11 "\t" $12 "\t" $13 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' | sed -e "s/_/:/g" | sed -e "s/BND/TRA/g" > ${TUMOUR}_${REFERENCE}.SV.bed
                            if [ `cat ${TUMOUR}_${REFERENCE}.SV.bed | wc -l` -gt 0 ]; then
                                    sed -i -e "1i BreakPoin1\tGene1\tBreakPoint2\tGene2\tType\tReadsRef\tReadsVar\tVAF" ${TUMOUR}_${REFERENCE}.SV.bed
                            else
                                    echo -e "BreakPoin1\tGene1\tBreakPoint2\tGene2\tType\tReadsRef\tReadsVar\tVAF" >> ${TUMOUR}_${REFERENCE}.SV.bed
                            fi
                            rm -r TEMP
                            mv ${TUMOUR}_${REFERENCE}.delly.vcf ${TUMOUR}_${REFERENCE}.SV.vcf
                            for j in `cut -f1 ${TUMOUR}_${REFERENCE}.SV.bed | sed -e "s/:/\t/" | tail -n+2`; do CHR=`echo ${j} | awk '{ print $1 }'`; BP=`echo ${j} | awk '{ print $2 }'`; TOTALREADS=`samtools view ${TUMOUR}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | wc -l`; XAPAREADS=`samtools view ${TUMOUR}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | grep XA: | grep -v pa: | grep "[;:]chr[0-9]*," | wc -l`; echo -e "${CHR}:${BP}\t${TOTALREADS}\t${XAPAREADS}"; done > XABP1
                            for j in `cut -f3 ${TUMOUR}_${REFERENCE}.SV.bed | sed -e "s/:/\t/" | tail -n+2`; do CHR=`echo ${j} | awk '{ print $1 }'`; BP=`echo ${j} | awk '{ print $2 }'`; TOTALREADS=`samtools view ${TUMOUR}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | wc -l`; XAPAREADS=`samtools view ${TUMOUR}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | grep XA: | grep -v pa: | grep "[;:]chr[0-9]*," | wc -l`; echo -e "${CHR}:${BP}\t${TOTALREADS}\t${XAPAREADS}"; done > XABP2
                            NVARS=`cat XABP1 | wc -l`
                            echo -e "Total\tXA\tXAPERC" > XA
                            for j in $( seq 1 $NVARS ); do XABP1N=`head -n ${j} XABP1 | tail -n1 | awk -F '\t' '{ print $2 }'`; XABP1XA=`head -n ${j} XABP1 | tail -n1 | awk -F '\t' '{ print $3 }'`; XABP2N=`head -n ${j} XABP2 | tail -n1 | awk -F '\t' '{ print $2 }'`; XABP2XA=`head -n ${j} XABP2 | tail -n1 | awk -F '\t' '{ print $3 }'`; let TOTAL=XABP1N+XABP2N; let TOTALXA=XABP1XA+XABP2XA; XAPERC=`echo -e "${TOTAL}\t${TOTALXA}" | awk '{ print $2 / $1 }'`; echo -e "${TOTAL}\t${TOTALXA}\t${XAPERC}"; done >> XA
                            paste ${TUMOUR}_${REFERENCE}.SV.bed XA > XAtemp && mv XAtemp ${TUMOUR}_${REFERENCE}.SV.bed
                            rm XA*
                    else
                            ## Somatic SV calling
                            echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Calling Tumour SVs" | tee -a ${LIST}.log
                            /data/Resources/Software/delly/src/delly call -x /data/Resources/Software/delly/excludeTemplates/human.hg38.excl.tsv -o ${TUMOUR}_${REFERENCE}.delly.bcf -g /data/Resources/References/${REFERENCE}/${REFERENCE}.fa ${TUMOUR}_${REFERENCE}.bam > /dev/null 2>&1
                            /data/Resources/Software/delly/src/delly merge -o sites.bcf ${TUMOUR}_${REFERENCE}.delly.bcf > /dev/null 2>&1
                            /data/Resources/Software/delly/src/delly call -g /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -v sites.bcf -o ${TUMOUR}_${REFERENCE}.delly.bcf -x /data/Resources/Software/delly/excludeTemplates/human.hg38.excl.tsv ${TUMOUR}_${REFERENCE}.bam > /dev/null 2>&1
                            bcftools view ${TUMOUR}_${REFERENCE}.delly.bcf | grep -v LowQual > ${TUMOUR}_${REFERENCE}.delly.vcf 2>&1
                            rm sites.bcf*
                            rm ${TUMOUR}_${REFERENCE}.delly.bcf*
                            mkdir TEMP
                            grep -v \# ${TUMOUR}_${REFERENCE}.delly.vcf | awk -F '\t' '{ print $1 "\t" $2 "\t" $8 "\t" $10}' | sed -e "s/\(\S\+\)\s\+\S\+SVTYPE=\(\w\+\)\S\+CHR2=\(\w\+\)\S\+END=\(\w\+\)\S\+/\1\t\3\t\4\t\2/" | sed -e "s/:/\t/g" | awk -F '\t' '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $14 + $16 "\t" $15 + $17 }' | awk -F '\t' '{ print $0 "\t" $7 / ($6 + $7) "\t" NR }' > TEMP/bed1
                            awk -F '\t' '{ print $1 "\t" $2 - 100 "\t" $2 + 100 "\t" $9 }' TEMP/bed1 | sort-bed - | bedmap --ec --echo --delim '\t' --echo-map-id --range 100 - /data/Resources/BEDs/All_Genes_hg38.unique.sorted2.broad.bed | awk -F '\t' '{ print $1 "_" $2 + 100 "\t" $5 "\t" $4 }' | sed -e "s/\t\t/\t\.\t/g" | sort -t$'\t' -k3,3g > TEMP/annot1
                            awk -F '\t' '{ print $3 "\t" $4 - 100 "\t" $4 + 100 "\t" $9 }' TEMP/bed1 | sort-bed - | bedmap --ec --echo --delim '\t' --echo-map-id --range 100 - /data/Resources/BEDs/All_Genes_hg38.unique.sorted2.broad.bed | awk -F '\t' '{ print $1 "_" $2 + 100 "\t" $5 "\t" $4 }' | sed -e "s/\t\t/\t\.\t/g" | sort -t$'\t' -k3,3g > TEMP/annot2
                            join -1 3 -2 3 TEMP/annot1 TEMP/annot2 | sed -e "s/\s\+/\t/g" | join -1 9 -2 1 TEMP/bed1 - | sed -e "s/\s\+/\t/g" | awk -F '\t' '{ print $10 "\t" $11 "\t" $12 "\t" $13 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' | sed -e "s/_/:/g" | sed -e "s/BND/TRA/g" > ${TUMOUR}_${REFERENCE}.SV.bed
                            if [ `cat ${TUMOUR}_${REFERENCE}.SV.bed | wc -l` -gt 0 ]; then
                                    sed -i -e "1i BreakPoin1\tGene1\tBreakPoint2\tGene2\tType\tReadsRef\tReadsVar\tVAF" ${TUMOUR}_${REFERENCE}.SV.bed
                            else
                                    echo -e "BreakPoin1\tGene1\tBreakPoint2\tGene2\tType\tReadsRef\tReadsVar\tVAF" >> ${TUMOUR}_${REFERENCE}.SV.bed
                            fi
                            rm -r TEMP
                            mv ${TUMOUR}_${REFERENCE}.delly.vcf ${TUMOUR}_${REFERENCE}.SV.vcf
                            for j in `cut -f1 ${TUMOUR}_${REFERENCE}.SV.bed | sed -e "s/:/\t/" | tail -n+2`; do CHR=`echo ${j} | awk '{ print $1 }'`; BP=`echo ${i} | awk '{ print $2 }'`; TOTALREADS=`samtools view ${TUMOUR}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | wc -l`; XAPAREADS=`samtools view ${TUMOUR}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | grep XA: | grep -v pa: | grep "[;:]chr[0-9]*," | wc -l`; echo -e "${CHR}:${BP}\t${TOTALREADS}\t${XAPAREADS}"; done > XABP1
                            for j in `cut -f3 ${TUMOUR}_${REFERENCE}.SV.bed | sed -e "s/:/\t/" | tail -n+2`; do CHR=`echo ${j} | awk '{ print $1 }'`; BP=`echo ${i} | awk '{ print $2 }'`; TOTALREADS=`samtools view ${TUMOUR}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | wc -l`; XAPAREADS=`samtools view ${TUMOUR}_${REFERENCE}.bam "${CHR}:${BP}-${BP}" 2> /dev/null | grep XA: | grep -v pa: | grep "[;:]chr[0-9]*," | wc -l`; echo -e "${CHR}:${BP}\t${TOTALREADS}\t${XAPAREADS}"; done > XABP2
                            NVARS=`cat XABP1 | wc -l`
                            echo -e "Total\tXA\tXAPERC" > XA
                            for j in $( seq 1 $NVARS ); do XABP1N=`head -n ${j} XABP1 | tail -n1 | awk -F '\t' '{ print $2 }'`; XABP1XA=`head -n ${j} XABP1 | tail -n1 | awk -F '\t' '{ print $3 }'`; XABP2N=`head -n ${j} XABP2 | tail -n1 | awk -F '\t' '{ print $2 }'`; XABP2XA=`head -n ${j} XABP2 | tail -n1 | awk -F '\t' '{ print $3 }'`; let TOTAL=XABP1N+XABP2N; let TOTALXA=XABP1XA+XABP2XA; XAPERC=`echo -e "${TOTAL}\t${TOTALXA}" | awk '{ print $2 / $1 }'`; echo -e "${TOTAL}\t${TOTALXA}\t${XAPERC}"; done >> XA
                            paste ${TUMOUR}_${REFERENCE}.SV.bed XA > XAtemp && mv XAtemp ${TUMOUR}_${REFERENCE}.SV.bed
                            rm XA*
                    fi
            done
        fi
        ## CNV calling
        if [ ${DATATYPE} != 'sWGS' ]; then
                if [ ${DATATYPE} = 'WGS' ] || [ ${DATATYPE} = 'sWGS' ]; then
                        BINLENGTH=1000
                else
                        BINLENGTH=0
                fi
                mkdir .cnv_temp
                ## preprocess intervals
                echo -e `date +%y%m%d\ %H%M%S\ %s` "     Preprocessing intervals" | tee -a ${LIST}.log
                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar PreprocessIntervals -L ${BED}.interval_list -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa --bin-length ${BINLENGTH} --interval-merging-rule OVERLAPPING_ONLY -O .cnv_temp/bed.interval_list > /dev/null 2>&1 
                ######################################################## NORMAL PON for READS
                ## Count Reads on Normals
                if [ `cat .AlignmentTemp_${LIST}/SamplesToCallGL.list | wc -l` -gt 0 ]; then
                        for i in `cat ${LIST}.pairs.txt | awk -F '\t' '{ print $1}'`; do
                                if [ -e ${i}_${REFERENCE}.bam ]; then
                                if [ ! -e .cnv_temp/${i}_${REFERENCE}.hdf5 ]; then
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Collecting read counts for PON" | tee -a ${LIST}.log
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CollectReadCounts -I ${i}_${REFERENCE}.bam -L .cnv_temp/bed.interval_list --interval-merging-rule OVERLAPPING_ONLY -O .cnv_temp/${i}_${REFERENCE}.hdf5 > /dev/null 2>&1 
                                        echo -e ".cnv_temp/${i}_${REFERENCE}.hdf5" >> .cnv_temp/PON.list
                                fi
                                fi
                        done
                        wait
                        ## Build the PON
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     Building the PON" | tee -a ${LIST}.log
                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CreateReadCountPanelOfNormals -I .cnv_temp/PON.list --minimum-interval-median-percentile 5.0 -O .cnv_temp/${LIST}.PON.hdf5 > /dev/null 2>&1
                else
                        BEDNOPATH=`echo $BED | sed -e "s/\/\S\+\///"`
                        if [ -e /data/Resources/PON/${BEDNOPATH}/PON.hdf5 ]; then
                                cp /data/Resources/PON/${BEDNOPATH}/PON.hdf5 .cnv_temp/${LIST}.PON.hdf5
                        else
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     No option available for the PON" | tee -a ${LIST}.log
                        fi
                fi
                ####################################################### END OF NORMAL PON for READS
                if [ -e .cnv_temp/${LIST}.PON.hdf5 ]; then
                        ## Count Reads on Tumours
                        for i in `cat ${LIST}.pairs.txt | awk -F '\t' '{ print $2}'`; do
                                if [ -e ${i}_${REFERENCE}.bam ]; then
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Collecting read counts of tumour sample" | tee -a ${LIST}.log
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CollectReadCounts -I ${i}_${REFERENCE}.bam -L .cnv_temp/bed.interval_list --interval-merging-rule OVERLAPPING_ONLY -O .cnv_temp/${i}_${REFERENCE}.hdf5 > /dev/null 2>&1 
                                fi
                        done
                        wait
                        ## Denoise sample data with PON and plot
                        for i in `cat ${LIST}.txt`; do
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Denoising sample using PON" | tee -a ${LIST}.log
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar DenoiseReadCounts -I .cnv_temp/${i}_${REFERENCE}.hdf5 --count-panel-of-normals .cnv_temp/${LIST}.PON.hdf5 --standardized-copy-ratios .cnv_temp/${i}_${REFERENCE}.standardizedCR.tsv --denoised-copy-ratios .cnv_temp/${i}_${REFERENCE}.denoisedCR.tsv > /dev/null 2>&1
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar PlotDenoisedCopyRatios --standardized-copy-ratios .cnv_temp/${i}_${REFERENCE}.standardizedCR.tsv --denoised-copy-ratios .cnv_temp/${i}_${REFERENCE}.denoisedCR.tsv --sequence-dictionary /data/Resources/References/${REFERENCE}/${REFERENCE}.dict --minimum-contig-length 46709983 --output .cnv_temp/plots --output-prefix ${i}_${REFERENCE} > /dev/null 2>&1
                        done
                        wait
                        ## Count alleles to use in combination with the cov to call CNVs in both somatic and germline
                        for i in `cat ${LIST}.pairs.txt`; do
                                ONGOING=`ps aux | grep CollectAllelicCounts | wc -l`
                                while [ ${ONGOING} -gt ${SIMSAM} ]; do
                                        ONGOING=`ps aux | grep CollectAllelicCounts | wc -l`
                                        echo -ne "\rProcessing ${ONGOING} samples ATM, waiting for some to finish"
                                        sleep 1
                                done
                                NORMAL=`echo ${i} | awk '{ print $1 }'`
                                TUMOUR=`echo ${i} | awk '{ print $2 }'`
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${TUMOUR}: Counting alleles in tumour sample" | tee -a ${LIST}.log
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CollectAllelicCounts -L .cnv_temp/bed.interval_list -I ${TUMOUR}_${REFERENCE}.bam -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -O .cnv_temp/${TUMOUR}_${REFERENCE}.allelicCounts.tsv > /dev/null 2>&1 &
                        done
                        wait
                        for i in `cat ${SOMATICFILE}.txt | awk -F '\t' '{ print $1}' | sort -u`; do
                                ONGOING=`ps aux | grep CollectAllelicCounts | wc -l`
                                while [ ${ONGOING} -gt ${SIMSAM} ]; do
                                        ONGOING=`ps aux | grep CollectAllelicCounts | wc -l`
                                        echo -ne "\rProcessing ${ONGOING} samples ATM, waiting for some to finish"
                                        sleep 1
                                done
                                if [ -e ${i}_${REFERENCE}.bam ]; then
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Counting alleles in normal sample" | tee -a ${LIST}.log
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CollectAllelicCounts -L .cnv_temp/bed.interval_list -I ${i}_${REFERENCE}.bam -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -O .cnv_temp/${i}_${REFERENCE}.allelicCounts.tsv > /dev/null 2>&1 &
                                fi
                        done
                        wait
                        ## Model segments using Cov and allele counts
                        for i in `cat ${LIST}.pairs.txt`; do
                                NORMAL=`echo ${i} | awk '{ print $1 }'`
                                TUMOUR=`echo ${i} | awk '{ print $2 }'`
                                if [ -e ${NORMAL}_${REFERENCE}.bam ]; then
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${TUMOUR}: Modelling segments in tumour sample" | tee -a ${LIST}.log
                                        java -Xmx100g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar ModelSegments --denoised-copy-ratios .cnv_temp/${TUMOUR}_${REFERENCE}.denoisedCR.tsv --allelic-counts .cnv_temp/${TUMOUR}_${REFERENCE}.allelicCounts.tsv --normal-allelic-counts .cnv_temp/${NORMAL}_${REFERENCE}.allelicCounts.tsv --output .cnv_temp/ --output-prefix ${TUMOUR}_${REFERENCE} > /dev/null 2>&1
                                else
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${TUMOUR}: Modelling segments in tumour sample" | tee -a ${LIST}.log
                                        java -Xmx100g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar ModelSegments --denoised-copy-ratios .cnv_temp/${TUMOUR}_${REFERENCE}.denoisedCR.tsv --allelic-counts .cnv_temp/${TUMOUR}_${REFERENCE}.allelicCounts.tsv --output .cnv_temp/ --output-prefix ${TUMOUR}_${REFERENCE} > /dev/null 2>&1
                                fi
                        done
                        wait
                        ## Call Segments
                        for i in `cat ${LIST}.pairs.txt`; do
                                NORMAL=`echo ${i} | awk '{ print $1 }'`
                                TUMOUR=`echo ${i} | awk '{ print $2 }'`
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${TUMOUR}: Calling CNVs in tumour sample" | tee -a ${LIST}.log
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CallCopyRatioSegments --input .cnv_temp/${TUMOUR}_${REFERENCE}.cr.seg --output .cnv_temp/${TUMOUR}_${REFERENCE}.called.seg > /dev/null 2>&1
                                grep -v @ .cnv_temp/${TUMOUR}_${REFERENCE}.called.seg | awk -F '\t' '{ if ( $5 >= 0.2 || $5 <= -0.2 ) print }' > .cnv_temp/${TUMOUR}_${REFERENCE}.CNV.bed
                        done
                        wait
                        ## Plot results
                        for i in `cat ${LIST}.pairs.txt`; do
                                NORMAL=`echo ${i} | awk '{ print $1 }'`
                                TUMOUR=`echo ${i} | awk '{ print $2 }'`
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${TUMOUR}: Plotting CNVs" | tee -a ${LIST}.log
                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar PlotModeledSegments --denoised-copy-ratios .cnv_temp/${TUMOUR}_${REFERENCE}.denoisedCR.tsv --allelic-counts .cnv_temp/${TUMOUR}_${REFERENCE}.hets.tsv -segments .cnv_temp/${TUMOUR}_${REFERENCE}.modelFinal.seg --sequence-dictionary /data/Resources/References/${REFERENCE}/${REFERENCE}.dict --minimum-contig-length 46709983 --output .cnv_temp/plots --output-prefix ${TUMOUR}_${REFERENCE} > /dev/null 2>&1
                                cp .cnv_temp/plots/${TUMOUR}_${REFERENCE}.modeled.png ${TUMOUR}_${REFERENCE}.CNV.png
                        done
                        wait
                        for i in `cat ${LIST}.pairs.txt`; do
                                NORMAL=`echo ${i} | awk '{ print $1 }'`
                                TUMOUR=`echo ${i} | awk '{ print $2 }'`
                                head -n1 .cnv_temp/${TUMOUR}_${REFERENCE}.CNV.bed > .cnv_temp/head
                                sed -i -e "s/\(.\+\)/\1\tGENES/" .cnv_temp/head
                                tail -n+2 .cnv_temp/${TUMOUR}_${REFERENCE}.CNV.bed | sort-bed - | bedmap --ec --echo --delim '\t' --echo-map-id --range 100 - /data/Resources/BEDs/All_Genes_hg38.unique.sorted2.broad.bed > .cnv_temp/tail
                                cat .cnv_temp/head .cnv_temp/tail > ${TUMOUR}_${REFERENCE}.CNV.bed
                        done
                        wait                
                else
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     No option available for CNV PON. Not calling somatic CNVs." | tee -a ${LIST}.log
                fi
                ## germline CNV calling for pairs
                touch .cnv_temp/germlinetocall.txt
                for i in `cat ${LIST}.pairs.txt | awk -F '\t' '{ print $1}' | sort -u`; do
                        if [ -e ${i}_${REFERENCE}.bam ]; then
                                echo -e "${i}" >> .cnv_temp/germlinetocall.txt
                        fi
                done
                if [ `cat .cnv_temp/germlinetocall.txt | wc -l` -gt 0 ]; then
                        ## preprocess intervals
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     Preprocessing intervals" | tee -a ${LIST}.log
                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar PreprocessIntervals -L ${BED}.interval_list -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa --bin-length ${BINLENGTH} --interval-merging-rule OVERLAPPING_ONLY -O .cnv_temp/bed.interval_list > /dev/null 2>&1 
                        ######################################################## NORMAL PON for READS
                        BEDNOPATH=`echo $BED | sed -e "s/\/\S\+\///"`
                        if [ -e /data/Resources/PON/${BEDNOPATH}/PON.hdf5 ]; then
                                cp /data/Resources/PON/${BEDNOPATH}/PON.hdf5 .cnv_temp/${LIST}.PON.hdf5
                        else
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     No option available for the PON" | tee -a ${LIST}.log
                        fi
                        ####################################################### END OF NORMAL PON for READS
                        if [ -e .cnv_temp/${LIST}.PON.hdf5 ]; then
                                ## Count Reads on Tumours
                                for i in `cat .cnv_temp/germlinetocall.txt`; do
                                        if [ -e ${i}_${REFERENCE}.bam ]; then
                                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Collecting read counts of normal sample" | tee -a ${LIST}.log
                                                java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CollectReadCounts -I ${i}_${REFERENCE}.bam -L .cnv_temp/bed.interval_list --interval-merging-rule OVERLAPPING_ONLY -O .cnv_temp/${i}_${REFERENCE}.hdf5 > /dev/null 2>&1 
                                        fi
                                done
                                wait
                                ## Denoise sample data with PON and plot
                                for i in `cat .cnv_temp/germlinetocall.txt`; do
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Denoising sample using PON" | tee -a ${LIST}.log
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar DenoiseReadCounts -I .cnv_temp/${i}_${REFERENCE}.hdf5 --count-panel-of-normals .cnv_temp/${LIST}.PON.hdf5 --standardized-copy-ratios .cnv_temp/${i}_${REFERENCE}.standardizedCR.tsv --denoised-copy-ratios .cnv_temp/${i}_${REFERENCE}.denoisedCR.tsv > /dev/null 2>&1
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar PlotDenoisedCopyRatios --standardized-copy-ratios .cnv_temp/${i}_${REFERENCE}.standardizedCR.tsv --denoised-copy-ratios .cnv_temp/${i}_${REFERENCE}.denoisedCR.tsv --sequence-dictionary /data/Resources/References/${REFERENCE}/${REFERENCE}.dict --minimum-contig-length 46709983 --output .cnv_temp/plots --output-prefix ${i}_${REFERENCE} > /dev/null 2>&1
                                done
                                wait
                                ## Count alleles to use in combination with the cov to call CNVs in both somatic and germline
                                for i in `cat .cnv_temp/germlinetocall.txt`; do
                                        ONGOING=`ps aux | grep CollectAllelicCounts | wc -l`
                                        while [ ${ONGOING} -gt ${SIMSAM} ]; do
                                                ONGOING=`ps aux | grep CollectAllelicCounts | wc -l`
                                                echo -ne "\rProcessing ${ONGOING} samples ATM, waiting for some to finish"
                                                sleep 1
                                        done
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Counting alleles in normal sample" | tee -a ${LIST}.log
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CollectAllelicCounts -L .cnv_temp/bed.interval_list -I ${i}_${REFERENCE}.bam -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa -O .cnv_temp/${i}_${REFERENCE}.allelicCounts.tsv > /dev/null 2>&1 &
                                done
                                wait
                                ## Model segments using Cov and allele counts
                                for i in `cat .cnv_temp/germlinetocall.txt`; do
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Modelling segments in normal sample" | tee -a ${LIST}.log
                                        java -Xmx100g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar ModelSegments --denoised-copy-ratios .cnv_temp/${i}_${REFERENCE}.denoisedCR.tsv --allelic-counts .cnv_temp/${i}_${REFERENCE}.allelicCounts.tsv --output .cnv_temp/ --output-prefix ${i}_${REFERENCE} > /dev/null 2>&1
                                done               
                                wait
                                ## Call Segments
                                for i in `cat .cnv_temp/germlinetocall.txt`; do
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Calling CNVs in normal sample" | tee -a ${LIST}.log
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar CallCopyRatioSegments --input .cnv_temp/${i}_${REFERENCE}.cr.seg --output .cnv_temp/${i}_${REFERENCE}.called.seg > /dev/null 2>&1
                                        grep -v @ .cnv_temp/${i}_${REFERENCE}.called.seg | awk -F '\t' '{ if ( $5 >= 0.2 || $5 <= -0.2 ) print }' > .cnv_temp/${i}_${REFERENCE}.CNV.bed
                                done
                                wait
                                ## Plot results
                                for i in `cat .cnv_temp/germlinetocall.txt`; do
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: Plotting CNVs" | tee -a ${LIST}.log
                                        java -Xmx20g -jar /data/Resources/Software/Javas/gatk-package-4.1.4.0-local.jar PlotModeledSegments --denoised-copy-ratios .cnv_temp/${i}_${REFERENCE}.denoisedCR.tsv   --allelic-counts .cnv_temp/${i}_${REFERENCE}.hets.tsv -segments .cnv_temp/${i}_${REFERENCE}.modelFinal.seg --sequence-dictionary /data/Resources/References/${REFERENCE}/${REFERENCE}.dict --minimum-contig-length 46709983 --output .cnv_temp/plots --output-prefix ${i}_${REFERENCE} > /dev/null 2>&1
                                        cp .cnv_temp/plots/${i}_${REFERENCE}.modeled.png ${i}_${REFERENCE}.CNV.png
                                done
                                wait
                                for i in `cat .cnv_temp/germlinetocall.txt`; do
                                        head -n1 .cnv_temp/${i}_${REFERENCE}.CNV.bed > .cnv_temp/head
                                        sed -i -e "s/\(.\+\)/\1\tGENES/" .cnv_temp/head
                                        tail -n+2 .cnv_temp/${i}_${REFERENCE}.CNV.bed | sort-bed - | bedmap --ec --echo --delim '\t' --echo-map-id --range 100 - /data/Resources/BEDs/All_Genes_hg38.unique.sorted2.broad.bed > .cnv_temp/tail
                                        cat .cnv_temp/head .cnv_temp/tail > ${i}_${REFERENCE}.CNV.bed
                                done
                                wait
                        else
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     No option available for CNV PON. Not calling germline CNVs." | tee -a ${LIST}.log
                        fi
                fi
                rm -r .cnv_temp/
        fi
        ####################################################### ichorCNA calls 
        if [ ${DATATYPE} = 'WGS' ] || [ ${DATATYPE} = 'sWGS' ]; then
                mkdir ichorCNA
                for i in `cat ${LIST}.pairs.txt | awk -F '\t' '{ print $2}'`; do
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: ichorCNA analysis" | tee -a ${LIST}.log
                        mkdir ichorCNA/${i}_${REFERENCE}/
                        cp ${i}_${REFERENCE}.bai ${i}_${REFERENCE}.bam.bai
                        /data/Resources/Software/hmmcopy_utils/bin/readCounter --window 1000000 --quality 20 --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" ${i}_${REFERENCE}.bam > ichorCNA/${i}_${REFERENCE}/${i}_${REFERENCE}.wig 2> /dev/null
                        Rscript /data/Resources/Software/ichorCNA/scripts/runIchorCNA.R --id ${i} --WIG ichorCNA/${i}_${REFERENCE}/${i}_${REFERENCE}.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 --gcWig /home/cebercoto/Downloads/ichorCNA/inst/extdata/gc_hg38_1000kb.wig --mapWig  /home/cebercoto/Downloads/ichorCNA/inst/extdata/map_hg38_1000kb.wig --centromere /home/cebercoto/Downloads/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt --normalPanel /home/cebercoto/Downloads/ichorCNA/inst/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" --estimateNormal True --estimatePloidy True --estimateScPrevalence True --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir ./ichorCNA/${i}_${REFERENCE}/ > /dev/null 2>&1
                        rm ${i}_${REFERENCE}.bam.bai 
                done
        else
                mkdir ichorCNA
                for i in `cat ${LIST}.pairs.txt | awk -F '\t' '{ print $2}'`; do
                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     ${i}: ichorCNA analysis" | tee -a ${LIST}.log
                        mkdir ichorCNA/${i}_${REFERENCE}/
                        cp ${i}_${REFERENCE}.bai ${i}_${REFERENCE}.bam.bai
                        /data/Resources/Software/hmmcopy_utils/bin/readCounter --window 1000000 --quality 20 --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" ${i}_${REFERENCE}.bam > ichorCNA/${i}_${REFERENCE}/${i}_${REFERENCE}.wig 2> /dev/null
                        Rscript /data/Resources/Software/ichorCNA/scripts/runIchorCNA.R --id ${i} --WIG ichorCNA/${i}_${REFERENCE}/${i}_${REFERENCE}.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 --gcWig /home/cebercoto/Downloads/ichorCNA/inst/extdata/gc_hg38_1000kb.wig --mapWig  /home/cebercoto/Downloads/ichorCNA/inst/extdata/map_hg38_1000kb.wig --centromere /home/cebercoto/Downloads/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt --normalPanel /home/cebercoto/Downloads/ichorCNA/inst/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" --estimateNormal True --estimatePloidy True --estimateScPrevalence True --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir ./ichorCNA/${i}_${REFERENCE}/ --exons.bed=${BED}.bed > /dev/null 2>&1
                        rm ${i}_${REFERENCE}.bam.bai 
                done
        fi
else
        echo -e `date +%y%m%d\ %H%M%S\ %s` "     --somatic ${SOMTAIC}: not recognized" | tee -a ${LIST}.log
fi
echo -e `date +%y%m%d\ %H%M%S\ %s` "     Finished calling variants" | tee -a ${LIST}.log