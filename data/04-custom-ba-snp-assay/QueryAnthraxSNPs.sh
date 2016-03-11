while read foldername
do
cat Final_Kmers_1793.txt | xargs -n 2 /data1/home/sandeep/Jelly_fish_2.0/bin/jellyfish query ${foldername}.jf > ${foldername}_Kmer_Counts.txt
done < list_2.txt
