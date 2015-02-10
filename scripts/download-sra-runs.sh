#! /bin/bash

mkdir sra-runs
cd sra-runs
prefetch -t ascp -a "/data1/home/rpetit/.aspera/connect/bin/ascp|/data1/home/rpetit/.aspera/connect/etc/asperaweb_id_dsa.openssh" --option-file ../data/SRP051511.txt 1> ../prefetch.out 2> ../prefetch.err
cd ..
