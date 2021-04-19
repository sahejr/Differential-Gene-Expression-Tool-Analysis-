#!/bin/bash

fastq-dump --split-files -O . SRX642055
seqtk sample -s 100 .fq 10000 > sub1.fq
seqtk sample -s 100 .fq 10000 > sub1.fq

fastq-dump --split-files -O . SRX642051
seqtk sample -s 100 read1.fq 10000 > sub1.fq
seqtk sample -s 100 read1.fq 10000 > sub1.fq
