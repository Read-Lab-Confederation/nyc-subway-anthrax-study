while read foldername
do
/data1/home/sandeep/Jelly_fish_2.0/bin/jellyfish count -m 31 -s 100M -t 12 -C ${foldername} -o ${foldername}.jf
done < list_2.txt 
