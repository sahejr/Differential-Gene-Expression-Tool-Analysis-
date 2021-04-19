#!/bin/bash

# Set variables used in each kallisto run
PROJECT=/home/linux/ieng6/cs185s/sdrandha/project
GTF=${PROJECT}/annotEdit.gtf
KINDEX=${PROJECT}/chr1Txn.idx

# Do a separate kallisto run for each dataset
for prefix in SRX642051 SRX642055 
do
    mkdir -p ${PROJECT}/${prefix}
    kallisto quant -t 3 -b 100 \
	-o ${PROJECT}/$prefix --gtf $GTF -i $KINDEX \
	${PROJECT}/reads/${prefix}_1.fq ${PROJECT}/reads/${prefix}_2.fq
done
