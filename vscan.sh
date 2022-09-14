#!/bin/bash

# automated virus discovery pipeline

#### CONFIGURATIONS ####
# 1. Software configs:
#   1.1 diamond executive
DIAMOND_EXEC="/usr/local/bin/diamond"
#   1.2 TransDecoder executive
TD_ORF_EXEC="TransDecoder.LongOrfs"
TD_PRED_EXEC="TransDecoder.Predict"
#   1.3 NCBI taxonomist executive
NCBI_NCBI_TAXONOMIST_EXEC="python3 /home/public/simbiont-js/tools/ncbi/ncbi.taxonomist.py"
#   1.4 csvtk, taxonkit, seqtk, seqkit executives
CSVTK_EXEC="/home1/panyuanfei/.conda/envs/yfpan/bin/csvtk"
TAXONKIT_EXEC="/home1/panyuanfei/.conda/envs/yfpan/bin/taxonkit"
SEQTK_EXEC="/home1/panyuanfei/.conda/envs/yfpan/bin/seqtk"
SEQKIT_EXEC="/home1/panyuanfei/.conda/envs/yfpan/bin/seqkit"
#   1.5 GNU Parallel executives
PARALLEL_EXEC="/usr/bin/parallel"

HMMSCAN_EXEC="/home1/panyuanfei/spack/opt/spack/linux-ubuntu18.04-zen2/aocc-3.2.0/hmmer-3.3.2-senbin2ucrtng4k6t76nwrkgoyxikvzy/bin/hmmscan"

# 2. Database configs
#   2.1 NCBI prot.accession2taxid database
PROT_ACCESSION2TAXID="/home/public/database/others/prot.accession2taxid"
#   2.2 diamond databases
DIAMOND_DB_NR="/shilab2/home/public/database/nr/nr.dmnd"
DIAMOND_DB_NR_VIRUS="/shilab2/home/panyuanfei/nr_virus_20210928.dmnd"
DIAMOND_DB_NR_RDRP="/shilab2/home/public/database/RdRp/mangRdRp20191008.dmnd"
DIAMOND_DB_VOG="/shilab2/home/panyuanfei/vog213.dmnd"
#   2.3 RVDB protein HMM profile
RVDB_PROT_HMM="/shilab2/home/panyuanfei/U-RVDBv24.1-prot.hmm"
#### CONFIGURATIONS ####

# user inputs
JOB_TITLE="18LSBat-R79"
CONTIG_MIN_LEN=1000
CONTIG_FASTA_FILE_RAW="${JOB_TITLE}.megahit.fa"
CONTIG_FASTA_FILE="${JOB_TITLE}.m${CONTIG_MIN_LEN}.megahit.fa"
MANIFEST_FILE="manifest.txt"
NUM_THREADS=40

# apply filter on to contig length
${SEQKIT_EXEC} seq -m ${CONTIG_MIN_LEN} ${CONTIG_FASTA_FILE_RAW} > ${CONTIG_FASTA_FILE}

# diamond blastx against nr
${DIAMOND_EXEC} blastx \
    -d ${DIAMOND_DB_NR_VIRUS} \
    -q ${CONTIG_FASTA_FILE} \
    -o ${JOB_TITLE}.blastx \
    -e 1E-3 \
    -k 1 \
    -p ${NUM_THREADS} \
    -f 6 qseqid qlen sseqid stitle pident length evalue \
    --more-sensitive

# get accession list from blastx result
cut -f3 ${JOB_TITLE}.blastx > ${JOB_TITLE}.acc

# tranlate accession into taxid
${PARALLEL_EXEC} --pipepart --block 400M \
    -j ${NUM_THREADS} \
    -a ${PROT_ACCESSION2TAXID} \
    grep -F -f ${JOB_TITLE}.acc > ${JOB_TITLE}.tax

# parse taxid
cat ${JOB_TITLE}.tax \
    | cut -f3 \
    | taxonkit lineage -c 2>/dev/null \
    | awk '$2>0' \
    | cut -f2- \
    | ${TAXONKIT_EXEC} reformat -F -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
    | cut -f1,3- \
    | ${CSVTK_EXEC} uniq -Ht \
    > ${JOB_TITLE}.tax.parsed

# merge parsed taxonomy with taxid table
${CSVTK_EXEC} join -Ht -L -f"3;1" ${JOB_TITLE}.tax ${JOB_TITLE}.tax.parsed \
    | awk '$3>0' \
    | cut -f"2,3,5-" \
    > ${JOB_TITLE}.tax.merged

# merge parsed taxonomy with blastx table
${CSVTK_EXEC} join -Ht -L -f"3;1" ${JOB_TITLE}.blastx ${JOB_TITLE}.tax.merged > ${JOB_TITLE}.blastx.tax

# remove previous results
if test ! -d "${JOB_TITLE}_COI"; then
    mkdir ${JOB_TITLE}_COI
else
    rm ${JOB_TITLE}_COI/*
fi
# select contigs of interest
IFS=$'\n'
for line in $(cat ${MANIFEST_FILE} | grep -v "^[[:space:]]*#" | sed '/^[[:space:]]*$/d'); do
    IFS=$'\t' read -r -a arr <<< "${line}"
    
    grp_name=${arr[0]}
    tax_level=${arr[1]}
    tax_names=${arr[2]}

    IFS=',' read -r -a tax_list <<< "${tax_names}"
    for tax in ${tax_list[@]}; do
        cat ${JOB_TITLE}.blastx.tax \
            | ${CSVTK_EXEC} grep -Ht -f ${tax_level} -p ${tax} \
            | cut -f1 \
            >> ${JOB_TITLE}_COI/${JOB_TITLE}.${grp_name}.list
        ${SEQTK_EXEC} subseq \
            ${CONTIG_FASTA_FILE} \
            ${JOB_TITLE}_COI/${JOB_TITLE}.${grp_name}.list \
            >> ${JOB_TITLE}_COI/${JOB_TITLE}.${grp_name}.fas
    done
done

# find all ORFs within contigs, then search against RdRp and nr_virus database
cd ${JOB_TITLE}_COI
for coi_fas in $(ls *.fas); do
    # TransDecoder find ORFs
    ${TD_ORF_EXEC} -t ${coi_fas} > /dev/null 2>&1
    ${TD_PRED_EXEC} -t ${coi_fas} > /dev/null 2>&1

    # blastp against RdRp database
    ${DIAMOND_EXEC} \
        blastp \
        -d ${DIAMOND_DB_NR_RDRP} \
        -q ${coi_fas}.transdecoder.pep \
        -o ${coi_fas}.rdrp.blastp \
        -e 1E-3 \
        -k 1 \
        -p ${NUM_THREADS} \
        -f 6 qseqid qlen sseqid stitle pident length evalue \
        --more-sensitive > /dev/null 2>&1

    # blastp against nr_virus database
    ${DIAMOND_EXEC} \
        blastp \
        -d ${DIAMOND_DB_NR_VIRUS} \
        -q ${coi_fas}.transdecoder.pep \
        -o ${coi_fas}.nr_virus.blastp.txt \
        -e 1E-3 \
        -p ${NUM_THREADS} \
        -f 6 qseqid qlen sseqid stitle pident length evalue \
        --more-sensitive > /dev/null 2>&1

    # blastp against VOG database
    ${DIAMOND_EXEC} \
        blastp \
        -d ${DIAMOND_DB_VOG} \
        -q ${coi_fas}.transdecoder.pep \
        -o ${coi_fas}.vog.blastp.txt \
        -e 1E-3 \
        -p ${NUM_THREADS} \
        -f 6 qseqid qlen sseqid stitle pident length evalue \
        --more-sensitive > /dev/null 2>&1

    # HMMER SCAN against RVDB
    ${HMMSCAN_EXEC} \
        -o ${coi_fas}.hmmer.out \
        --tblout ${coi_fas}.hmmer.tbl \
        --domtblout ${coi_fas}.hmmer.domtbl \
        --pfamtblout ${coi_fas}.hmmer.pfamtbl \
        ${RVDB_PROT_HMM} \
        ${coi_fas}.transdecoder.pep
done
cd ..