#!/usr/bin/env python
import re,string,time
import numpy as np
import sys

if len(sys.argv) != 7:
	print('Usage: ./contact_processing.py contacts_file output_file num_frames'+'\n'+
				'                               num_residues threshold num_chains'+ '\n\n\t' +  
									 '- num_residues is num amino acids in chain' + '\n\t' + 
									 '- threshold is percent of max contact to use to define strong' + '\n\t' +
									 '  contacts'+'\n\n\t' + 
									 '  The output file lists for each residue the percentage of time it was'+'\n\t'
									 '  in strong contact, as defined by the threshold')
	exit()

# contact file output from vmd
vmd_contacts_file=sys.argv[1]
# output file of average normalized contact per residue
output_file=sys.argv[2]
# number of frames in trajectory used in vmd contact tool
numframes=int(sys.argv[3])+2
# number of residues on a chain
numresidues=int(sys.argv[4])
# threshold percentage of max contact, above which we consider a strong contact
threshold=float(sys.argv[5])
# number of chains
numunits=int(sys.argv[6])

# open both files
with open(vmd_contacts_file,'r') as file1, open(output_file,'w') as file2:
	# split with regexp the file by line
	entries = re.split('\n+', file1.read())
	# get length of data per chain over traj
	total=numframes*numresidues  
	# initialize container for contact data for each residue, 
  # averaged over chains and frames
	allData = np.empty((0,numframes-2), int) 
	# loop over entries, skipping every numframes
	# i starts at 9 to skip header in vmd_contacts_file
	for i in range(9,len(entries),numframes):
		# check if we're on a line for particular res 
		try:
			go=False
			for j in range(0,numunits-1):
				if entries[i+j*total]==entries[i+(j+1)*total]:
					go=True
				else:
					go=False
					break
			# if we are, get data for that res accross chains and frames
			if go:
				A = np.zeros(numframes-2, int)	
				for j in range(0,numunits-1):
					A = A + np.array([int(g[-2:]) for g in entries[i+2+j*total:i+(j*total)+numframes]])
				# populate allData with average data for the res
				allData = np.concatenate((allData,[A/float(numunits)]))
			i=i+numframes
		except IndexError:
			break

	# compare the contact data with a % of the max and 
	# sum up a binary value to count how many frames a res was in strong contact
	# the final result is percentage of trajectory in which a residue is in strong
	# contact
	for i in allData:
		compare = i >= threshold*max(i)
		file2.write(str(sum(compare.astype(int))/float(len(i)))+'\n')
file1.close()
file2.close()
