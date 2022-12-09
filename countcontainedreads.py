#Usage : python <file> <sorted-bed-file>
#Sorting must be done using first column as key in ascending order, followed by second column as key in ascending order, followed by third column as key in descending order

import sys
import csv

containedReads = 0
containedReadsAboveCutoffLen = 0
totalReads = 0
currentChromosome = ""
rightMostOffset = 0
minimumLengthCutoff = 0

if len(sys.argv) >= 3:
    minimumLengthCutoff = int(sys.argv[2])

with open(sys.argv[1], 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        startOffset = int(row[1])
        endOffset = int(row[2])
        if (row[0] != currentChromosome):
            currentChromosome = row[0]
            rightMostOffset = endOffset
        else:
            if endOffset > rightMostOffset:
                rightMostOffset = endOffset
            else:
                containedReads =  containedReads + 1
                if endOffset - startOffset >= minimumLengthCutoff:
                    containedReadsAboveCutoffLen = containedReadsAboveCutoffLen + 1;

        totalReads = totalReads + 1

print("Total reads = ", totalReads)
print("Contained reads = ", containedReads)
print("Contained reads = ", containedReadsAboveCutoffLen, " with length >=", minimumLengthCutoff)
