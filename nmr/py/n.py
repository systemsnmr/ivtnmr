# Open/create file with notes to current NMR dataset
# Needs to be run from an *experiment* window in TopSpin (i.e. will not work if called from edte/stdisp/etc)
# @Yaroslav 2012

# import routines to run shell commands
from subprocess import call
import os.path

# get info about current dataset
curdat = CURDATA()

# create path to notes file, using dataset name and TopSpin data folder
datasetNotes = curdat[3]+"/"+curdat[0]+"/notes.txt"

# If the file doesn't exist yet - create one, appending dataset title to the top.
if not(os.path.isfile(datasetNotes)):
	datasetHeader = "\n## "+curdat[0]+"\n===================================="
	f = open(datasetNotes, 'w')
	f.write( datasetHeader )
	f.close()

# open text editor with notes file
call(["gedit", datasetNotes])

#MSG(datasetNotes)
#MSG(curdat[0]+", "+ curdat[1]+", "+curdat[2]+", "+curdat[3])
#MSG("name="+curdat[0] + ", procno=" +curdat[2])