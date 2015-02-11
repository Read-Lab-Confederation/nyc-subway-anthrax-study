#! /bin/bash

while read line
do

prefetch -a "/data1/home/rpetit/.aspera/connect/bin/ascp|/data1/home/rpetit/.aspera/connect/etc/asperaweb_id_dsa.openssh" $line

done < $1
