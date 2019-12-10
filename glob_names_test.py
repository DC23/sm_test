# Create a list of files for use in the pipeline. Also needs to pair the sample names up and run them into paired end reads.

import glob
import os

raw_file_loc = '/OSM/CBR/AF_DATASCHOOL/input/2019-04-12_Transcritome'

raw_file_names = [os.path.basename(x) for x in glob.glob(raw_file_loc + '/*1.fastq.gz')[0:-7]]
#raw_file_names = raw_file_names[:-7]
print(raw_file_names)
print("Done")

munged_file_names = []
for i in raw_file_names:
    munged_file_names.append(i[0:-12])

print(munged_file_names)