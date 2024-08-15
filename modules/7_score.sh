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
CORES2=`echo "${CORES} * ${SIMSAM}" | bc -l`
IFS=$'\n'
if [ ${SOMATIC} = 'false' ]; then
        if [ -e ${BAMPATH}.bampath ]; then
                for i in `cat ${LIST}.txt`; do
                        ONGOING=`ps aux | grep BASS | wc -l`
                        while [ ${ONGOING} -gt ${CORES2} ]; do
                                ONGOING=`ps aux | grep BASS | wc -l`
                                echo -ne "\r${ONGOING} scores going on ATM, waiting for some to finish"
                                sleep 1
                        done
                        THISBAMPATH=`grep -w ${i} ${BAMPATH}.bampath | awk '{ print $2 }'`
                        if [ ! -e ${THISBAMPATH}/${i}_${REFERENCE}.bam ]; then
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     WARNING: sample ${i} bam file does not exist, not correcting for multimaps" | tee -a ${LIST}.log
                        fi
                        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Annotating and scoring individual ${i}" | tee -a ${LIST}.log
                        BASS --target ${i}_${REFERENCE}.QC -M individual -m ${BASSMODEL} -b ${REFERENCE} --MAF ${BASSMAF} --bam ${THISBAMPATH}/${i}_${REFERENCE} >/dev/null 2>&1 &
                done
                wait
        else
                for i in `cat ${LIST}.txt`; do
                        ONGOING=`ps aux | grep BASS | wc -l`
                        while [ ${ONGOING} -gt ${CORES2} ]; do
                                ONGOING=`ps aux | grep BASS | wc -l`
                                echo -ne "\r${ONGOING} scores going on ATM, waiting for some to finish"
                                sleep 1
                        done
                        if [ ! -e ${i}_${REFERENCE}.bam ]; then
                                echo -e `date +%y%m%d\ %H%M%S\ %s` "     WARNING: sample ${i} bam file does not exist, not correcting for multimaps" | tee -a ${LIST}.log
                        fi
                        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Annotating and scoring individual ${i}" | tee -a ${LIST}.log
                        BASS --target ${i}_${REFERENCE}.QC -M individual -m ${BASSMODEL} -b ${REFERENCE} --MAF ${BASSMAF} --bam ${i}_${REFERENCE} >/dev/null 2>&1 &
                done
                wait
        fi
        wait
        for i in `cat ${LIST}.txt`; do
                ( head -n1 ${i}_${REFERENCE}.QC.scores.txt & tail -n+2 ${i}_${REFERENCE}.QC.scores.txt | sed -e "s/\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\(.\+\)/\1_\2-\3\t\1\t\2\t\3\4/" | sort -t$'\t' -k1,1 -k19,19nr | sort -t$'\t' -u -k1,1 | sort -t$'\t' -k32,32nr | cut -f1 --complement ) > ${i}_${REFERENCE}.QC.scores.filter.txt              
        done
        rm *.QC.log
        rm *.hom
elif [ ${SOMATIC} = 'true' ]; then
        if [ -e ${BAMPATH}.bampath ]; then
                for i in `cat ${SOMATICFILE}.txt`; do
                        ONGOING=`ps aux | grep BASS | wc -l`
                        while [ ${ONGOING} -gt ${CORES} ]; do
                                ONGOING=`ps aux | grep BASS | wc -l`
                                echo -ne "\r${ONGOING} scores going on ATM, waiting for some to finish"
                                sleep 1
                        done
                        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Annotating and scoring individual ${i}" | tee -a ${LIST}.log
                        NORMAL=`echo ${i} | awk '{ print $1 }'`
                        TUMOUR=`echo ${i} | awk '{ print $2 }'`
                        if [ -e ${NORMAL}_${REFERENCE}.bam ]; then
                                THISBAMPATH=`grep -w ${TUMOUR} ${BAMPATH}.bampath | awk '{ print $2 }'`
                                if [ ! -e ${THISBAMPATH}/${TUMOUR}_${REFERENCE}.bam ]; then
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     WARNING: sample ${TUMOUR} bam file does not exist, not correcting for multimaps" | tee -a ${LIST}.log
                                fi
                                BASS --target ${TUMOUR}_${REFERENCE}_somatic.QC --somatic true -b ${REFERENCE} --MAF ${BASSMAF} --bam ${THISBAMPATH}/${TUMOUR}_${REFERENCE} >/dev/null 2>&1 &
                                BASS --target ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC --somatic true -b ${REFERENCE} --MAF ${BASSMAF} --bam ${THISBAMPATH}/${TUMOUR}_${REFERENCE} >/dev/null 2>&1 &
                                if [ ! -e .Score_${NORMAL}_${REFERENCE}.QC_TempFiles ]; then
                                        BASS --target ${NORMAL}_${REFERENCE}.QC -M individual -m ${BASSMODEL} -b ${REFERENCE} --MAF ${BASSMAF} --bam ${THISBAMPATH}/${NORMAL}_${REFERENCE} >/dev/null 2>&1 &
                                fi
                        else
                                THISBAMPATH=`grep -w ${TUMOUR} ${BAMPATH}.bampath | awk '{ print $2 }'`
                                if [ ! -e ${THISBAMPATH}/${TUMOUR}_${REFERENCE}.bam ]; then
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     WARNING: sample ${TUMOUR} bam file does not exist, not correcting for multimaps" | tee -a ${LIST}.log
                                fi
                                BASS --target ${TUMOUR}_${REFERENCE}_somatic.QC --somatic true -b ${REFERENCE} --MAF ${BASSMAF} --bam ${THISBAMPATH}/${TUMOUR}_${REFERENCE} >/dev/null 2>&1 &
                                BASS --target ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC --somatic true -b ${REFERENCE} --MAF ${BASSMAF} --bam ${THISBAMPATH}/${TUMOUR}_${REFERENCE} >/dev/null 2>&1 &
                        fi
                done
                wait
        else
                for i in `cat ${SOMATICFILE}.txt`; do
                        ONGOING=`ps aux | grep BASS | wc -l`
                        while [ ${ONGOING} -gt ${CORES} ]; do
                                ONGOING=`ps aux | grep BASS | wc -l`
                                echo -ne "\r${ONGOING} scores going on ATM, waiting for some to finish"
                                sleep 1
                        done
                        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     Annotating and scoring individual ${i}" | tee -a ${LIST}.log
                        NORMAL=`echo ${i} | awk '{ print $1 }'`
                        TUMOUR=`echo ${i} | awk '{ print $2 }'`
                        if [ -e ${NORMAL}_${REFERENCE}.bam ]; then
                                if [ ! -e ${TUMOUR}_${REFERENCE}.bam ]; then
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     WARNING: sample ${TUMOUR} bam file does not exist, not correcting for multimaps" | tee -a ${LIST}.log
                                fi
                                BASS --target ${TUMOUR}_${REFERENCE}_somatic.QC --somatic true -b ${REFERENCE} --MAF ${BASSMAF} --bam ${TUMOUR}_${REFERENCE} >/dev/null 2>&1 &
                                BASS --target ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC --somatic true -b ${REFERENCE} --MAF ${BASSMAF} --bam ${TUMOUR}_${REFERENCE} >/dev/null 2>&1 &
                                if [ ! -e .Score_${NORMAL}_${REFERENCE}.QC_TempFiles ]; then
                                        BASS --target ${NORMAL}_${REFERENCE}.QC -M individual -m ${BASSMODEL} -b ${REFERENCE} --MAF ${BASSMAF} --bam ${NORMAL}_${REFERENCE} >/dev/null 2>&1 &
                                fi
                        else
                                if [ ! -e ${TUMOUR}_${REFERENCE}.bam ]; then
                                        echo -e `date +%y%m%d\ %H%M%S\ %s` "     WARNING: sample ${TUMOUR} bam file does not exist, not correcting for multimaps" | tee -a ${LIST}.log
                                fi
                                BASS --target ${TUMOUR}_${REFERENCE}_somatic.QC --somatic true -b ${REFERENCE} --MAF ${BASSMAF} --bam ${TUMOUR}_${REFERENCE} >/dev/null 2>&1 &
                                BASS --target ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC --somatic true -b ${REFERENCE} --MAF ${BASSMAF} --bam ${TUMOUR}_${REFERENCE} >/dev/null 2>&1 &
                        fi
                done
                wait
        fi
        ## filtered annotations: apply MuTect2 filters and keep only the most cannical transcript affected by the mutation
        for i in `cat ${SOMATICFILE}.txt`; do
                NORMAL=`echo ${i} | awk '{ print $1 }'`
                TUMOUR=`echo ${i} | awk '{ print $2 }'`
                if [ -e ${NORMAL}_${REFERENCE}.bam ]; then
                        ( head -n1 ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.scores.txt & tail -n+2 ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.scores.txt | sed -e "s/\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\(.\+\)/\1_\2-\3\t\1\t\2\t\3\4/" | sort -t$'\t' -k1,1 -k19,19nr | sort -t$'\t' -u -k1,1 | sort -t$'\t' -k33,33nr | cut -f1 --complement | awk -F '\t' '{ if ( $31 == "-9" || $31 == "clustered_events" || $31 == "t_lod" || $31 == "clustered_events;t_lod" ) print }' | awk -F '\t' '{ if ( $30 > 200 ) print }' | awk -F '\t' '{ if ( $29 > 0.01 ) print }' ) > ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.scores.filter.txt
                        ( head -n1 ${TUMOUR}_${REFERENCE}_somatic.QC.scores.txt & tail -n+2 ${TUMOUR}_${REFERENCE}_somatic.QC.scores.txt | sed -e "s/\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\(.\+\)/\1_\2-\3\t\1\t\2\t\3\4/" | sort -t$'\t' -k1,1 -k19,19nr | sort -t$'\t' -u -k1,1 | sort -t$'\t' -k33,33nr | cut -f1 --complement | awk -F '\t' '{ if ( $31 == "-9" || $31 == "clustered_events" || $31 == "t_lod" || $31 == "clustered_events;t_lod" ) print }' | awk -F '\t' '{ if ( $30 > 200 ) print }' | awk -F '\t' '{ if ( $29 > 0.01 ) print }' ) > ${TUMOUR}_${REFERENCE}_somatic.QC.scores.filter.txt
                        ( head -n1 ${NORMAL}_${REFERENCE}.QC.scores.txt & tail -n+2 ${NORMAL}_${REFERENCE}.QC.scores.txt | sed -e "s/\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\(.\+\)/\1_\2-\3\t\1\t\2\t\3\4/" | sort -t$'\t' -k1,1 -k19,19nr | sort -t$'\t' -u -k1,1 | sort -t$'\t' -k32,32nr | cut -f1 --complement ) > ${NORMAL}_${REFERENCE}.QC.scores.filter.txt
                else
                        ( head -n1 ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.scores.txt & tail -n+2 ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.scores.txt | sed -e "s/\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\(.\+\)/\1_\2-\3\t\1\t\2\t\3\4/" | sort -t$'\t' -k1,1 -k19,19nr | sort -t$'\t' -u -k1,1 | sort -t$'\t' -k33,33nr | cut -f1 --complement | awk -F '\t' '{ if ( $31 == "-9" || $31 == "clustered_events" || $31 == "t_lod" || $31 == "clustered_events;t_lod" ) print }' | awk -F '\t' '{ if ( $30 > 200 ) print }' | awk -F '\t' '{ if ( $29 > 0.01 ) print }' ) > ${TUMOUR}_${REFERENCE}_somaticNOGLF.QC.scores.filter.txt
                        ( head -n1 ${TUMOUR}_${REFERENCE}_somatic.QC.scores.txt & tail -n+2 ${TUMOUR}_${REFERENCE}_somatic.QC.scores.txt | sed -e "s/\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\(.\+\)/\1_\2-\3\t\1\t\2\t\3\4/" | sort -t$'\t' -k1,1 -k19,19nr | sort -t$'\t' -u -k1,1 | sort -t$'\t' -k33,33nr | cut -f1 --complement | awk -F '\t' '{ if ( $31 == "-9" || $31 == "clustered_events" || $31 == "t_lod" || $31 == "clustered_events;t_lod" ) print }' | awk -F '\t' '{ if ( $30 > 200 ) print }' | awk -F '\t' '{ if ( $29 > 0.01 ) print }' ) > ${TUMOUR}_${REFERENCE}_somatic.QC.scores.filter.txt
                fi
        done
        ## merge duplicate annotations
        if [ ${DUPLICATES} = 'true' ]; then
                for i in `cat ${DUPFILE}.txt`; do
                        if [ `echo ${i} | awk -F'\t' '{print NF; exit}'` == 2 ]; then
                                echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     merging duplicates ${i}" | tee -a ${LIST}.log 
                                DUP1=`echo ${i} | awk '{ print $1 }'`
                                DUP2=`echo ${i} | awk '{ print $2 }'`
                                head -n1 ${DUP1}_${REFERENCE}_somatic.QC.scores.filter.txt > .DUPLICATES.dup1.head
                                tail -n+2 ${DUP2}_${REFERENCE}_somatic.QC.scores.filter.txt | awk -F '\t' '{ print $1 "_" $2 "-" $3 }' > .DUPLICATES.dup2.variants
                                tail -n+2 ${DUP1}_${REFERENCE}_somatic.QC.scores.filter.txt | awk -F '\t' '{ print $1 "_" $2 "-" $3 "\t" $0 }' > .DUPLICATES.dup1.tofilter
                                (cat .DUPLICATES.dup1.head & grep -w -f .DUPLICATES.dup2.variants .DUPLICATES.dup1.tofilter ) > ${DUP1}-${DUP2}_${REFERENCE}_somatic.QC.scores.filter.txt
                        else
                                echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     ${i}: exactly 2 duplicates not provided, skipping" | tee -a ${LIST}.log 
                        fi
                done
        fi
        rm .log
        rm *.QC.log
else
        echo -e "\r"`date +%y%m%d\ %H%M%S\ %s` "     --somatic ${SOMATIC}: option not recognized" | tee -a ${LIST}.log
        exit
fi
wait
