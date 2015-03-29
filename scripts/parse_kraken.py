#!/usr/bin/python

import sys
input1 = sys.argv[1]
input2 = sys.argv[2]
count = 0
anthracis = 0
bacillus = 0
bacteria = 0

fh1 = open(input1, 'r')
fh2 = open(input2, 'w')
#print >> fh2, "Bacillus_anthracis" + "\t" + "Other_Bacillus_sp" + "\t" + "Other_bacteria"
for line in fh1:
    fields = line.rstrip( "\r\n" ).split( "\t" )
    count += 1
    if count > 10:
       if fields[3] == 'S' or '-':
          if fields[5].find("Bacillus anthracis") != -1:
             anthracis += int(fields[1])
          elif fields[5].find("Bacillus") != -1:
              bacillus += int(fields[1])
          else :
              bacteria += int(fields[1])
print >> fh2, str(anthracis) + "\n" + str(bacillus) + "\n" + str(bacteria)
fh1.close()
fh2.close()