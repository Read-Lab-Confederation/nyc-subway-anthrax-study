while read foldername
do
mkdir ${foldername}
cd ${foldername}
/data1/home/sandeep/KRAKEN/kraken --db /data1/home/sandeep/KRAKEN/BACTERIA --gzip-compressed /data1/home/rpetit/readgp/nyc-subway-metagenome/sra-pathogens/yersinia/${foldername}/*_1.fastq.gz /data1/home/rpetit/readgp/nyc-subway-metagenome/sra-pathogens/yersinia/${foldername}/*_2.fastq.gz --paired --fastq-input --threads 20 --classified-out ${foldername}_Classified_Reads.fastq > KRACKEN_${foldername}.txt
/data1/home/sandeep/KRAKEN/kraken-report --db /data1/home/sandeep/KRAKEN/BACTERIA KRACKEN_${foldername}.txt > FINAL_KRACKEN_REPORT_${foldername}.txt
cd ..
done < /data1/home/sandeep/KRACKEN_NYC/YERSINIA/list_3.txt  
