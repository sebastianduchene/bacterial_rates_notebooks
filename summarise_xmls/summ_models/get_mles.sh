#!/bin/bash

log_files=`ls | awk '/mle[.]log/'`

for i in $log_files; do
    echo wrking on $i
    cp $i temp_mle.log
    echo $i >> out_data.txt
#    ~/BEASTv1.8.1/bin/beast temp_mle.xml > temp_result.txt
	~/Downloads/beast1/bin/beast temp_mle.xml > temp_result.txt
    awk '/pathLikelihood/' temp_result.txt >> out_data.txt
    rm temp_result.txt
    rm temp_mle.log
done

