import os
from shutil import copytree, ignore_patterns

raw_direc = 'DATA'
filtered_direc = 'DATA_FILTERED'

#Copies over every *.jpg (spin evolution profile, avalanche distribution plot, comparative avalanche distribution plot, and any other misc images) and other avalanche distribution .txt files (which don't fall into the ignorable pattern) to the corresponding directory in 'DATA, FILTERED'
#Basically, this should include every file except the raw spin data files (the vast majority of the used storage), any python scripts, and pycache
working_direc = os.getcwd()
copytree(working_direc + '/' + raw_direc, working_direc + '/' + filtered_direc, ignore=ignore_patterns('?.txt', '??.txt', '???.txt', '????.txt', '*MemRight.txt', '*MemUp.txt', '*.py', '__*')) #can add more '??...?.txt' if you're got > 9999 instances per ensemble

#Looks through each subdirectory in and finds the total number of instances considered in this ensemble
#Creates a number.txt file with the number of instances written in the file (so it's clear how many instances were simulated, without needing the raw data)
for subdirec in os.listdir(working_direc + '/' + raw_direc):
    num_of_instances = 0
    if subdirec[0:2] == 'sz': #only considers directories of a particular size (others have no '?.txt' files)
        for file in os.listdir(working_direc + '/' + raw_direc + '/' + subdirec):
            if len(file[:-4]) <= 4 and file[-4:] == '.txt': #filters for '?.txt', '??.txt', '???.txt', or '????.txt' files
                num_of_instances += 1
    ###print(subdirec + ' contains ' + str(num_of_instances) + ' instances')

        instance_num_file = open(working_direc + '/' + filtered_direc + '/' + subdirec + '/instance_num.txt', 'w')
        instance_num_file.write(str(num_of_instances))

print('Directories copied!')
