# Automatically read/overlay a range of spectra in the .md mode
# Script prompts user for the first and last exp in the to be read.
# If providing a RANGE of expts - just give first and last - WITHOUT A DASH!

# Automatically skips empty experiments  w/o complaint.
# @Yaroslav Nikolaev. 2013-2017

import os

# for checking if expts exist
curdat = CURDATA() # get current dataset
datasetFolder = curdat[3]+"/"+curdat[0]+"/"

####### Actual script #######
if len(sys.argv) == 3:
	first = int(sys.argv[1])
	last = int(sys.argv[2])
elif len(sys.argv) == 4:
	first = int(sys.argv[1])
	last = int(sys.argv[2])
	increment = int(sys.argv[3])
else:
	result = INPUT_DIALOG("", "Please enter first and last expnos to overlay:",
	["first = ", "last = "], ["1", "1"])

	if result == None:
		EXIT()
	
	first = int(result[0])
	last = int(result[1])

if len(sys.argv) == 4:
	for expno in range(first,last+1,increment): # Loop to copy every experiment in the block
		# Read expt only if it exists (skips empty expts):
		if os.path.exists(datasetFolder+"/"+str(expno)):
			XCMD('re '+str(expno), wait = WAIT_TILL_DONE) # read
		
else:
	for expno in range(first,last+1): # Loop to copy every experiment in the block
		# Read expt only if it exists (skips empty expts):
		if os.path.exists(datasetFolder+"/"+str(expno)):
			XCMD('re '+str(expno), wait = WAIT_TILL_DONE) # read