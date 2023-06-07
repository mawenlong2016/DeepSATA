#!/bin/bash

if [ "${1}" == "train" ];then  
	grep 'MOTIF' ${2} | awk '{print $2}' | uniq > motif.name
	if [ ! -d ./${1}_motif ]
	then mkdir -p ./${1}_motif
	fi
	Rscript ../../code/7_filterTF.R ${3} ${1}/fimo_back.tsv	${1}
	cat ${1}_motif/*.count > motif.result
	sort -n -r -k 1 motif.result | head -10 | awk '{print $2}' > motif.enrich.name
	Rscript ../../code/7_filterFimo.R ${3} ${1}/fimo_back.tsv ${1}	
	cat ${1}_motif/*.fimo > ${1}/fimo.tsv
	rm -r ${1}_motif	
elif [ "${1}" == "valid" ];then
	if [ ! -d ./${1}_motif ]
	then mkdir -p ./${1}_motif
	fi
	Rscript ../../code/7_filterFimo.R ${3} ${1}/fimo_back.tsv ${1}	
	cat ${1}_motif/*.fimo > ${1}/fimo.tsv
	rm -r ${1}_motif 
elif [ "${1}" == "test" ];then
	if [ ! -d ./${1}_motif ]
	then mkdir -p ./${1}_motif
	fi
	Rscript ../../code/7_filterFimo.R ${3} ${1}/fimo_back.tsv ${1}	
	cat ${1}_motif/*.fimo > ${1}/fimo.tsv
	rm -r ${1}_motif
else
    if [ ! -d ./${1}_motif ]
    then mkdir -p ./${1}_motif
    fi
    Rscript ../../code/7_filterFimo.R ${3} ${1}/fimo_back.tsv ${1}
    cat ${1}_motif/*.fimo > ${1}/fimo.tsv
    rm -r ${1}_motif
fi



