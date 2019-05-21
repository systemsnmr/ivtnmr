# Automatically create series of repeating datasets for time-resolved NMR.
# For example if you have [1D-31P, tocsy, 1D-1H] and want to run this set of experiments for 5 hours, 
# but do not want to manually do 'edc' few dozen times.
# Required settings/parameters:
# 1) Reference set of experiments to be repeated (time=0): list of exp numbers in certain order.
# 2) Set of experiment numbers for the time-point #1: list of experiments (will duplicate Reference set & blank the expt titles). 
# 3) (rough) times of individual experiments
# 4) (by default this is now calculated automatically): experiment number at which to start the series
# 5) time period to run the series
# To be executed from within a dataset (not EDTE, etc).
# @Yaroslav Nikolaev. 2013-..
# v002 (2015-06-23): included flag/option to also copy initial set of expts + cleanup their titles
# (2018-06-20): modified to use "setup_case" switching, instead of many if-elses.

#### TODO ####
# - functionalize the script (main(), etc)
# - Include a warning-check to tell in which dataset the series are going to be created! (often doing in wrong set :)
from TopCmds import *
from subprocess import call
#import sys, os, time, shutil, commands, collections
import os, time
import datetime
import sys # to get own filename
import shutil # to copy files

### SETTINGS
#===============================
setup_case_list = [
	'HN_31P_iminos_1H_tocsy_HN_31P_iminos_1H', #0		# IVTNMR -- WITH HSQC, tocsy every second round
	'31P_iminos_1H_tocsy_31P_iminos_1H', 			 #1		# IVTNMR -- no HSQC, tocsy every second round
	'mapper_many_hsqcs', 										   #2		# mapper -- MANY_HSQCs 
	'mapper_many_hsqcs_and_tocsy', 						 #3		# mapper -- MANY_HSQCs + TOCSY
	'SNIC_no_15Nprot', 												 #4
	'SNIC_with_15Nprot', 											 #5
	'5Pz', 																		 #6   # HSQCs a bit more frequent than other
	'HN_31P_every_second_1D1H',							   #7   # HN-P-H-HN-P (for 20-30uM protein)
	'HN_31P_every_second_1D1H_or_sofast',		   #8   # HN-P-H-HN-P-SF (for 50uM protein)
	]
	
setup_case = setup_case_list[2] # INDEXING STARTS FROM 0!! select item from above list

timeToRun = 24*60 # defined in minutes

# Transcription experiments can be recorded in two modes:
# CO-transcriptional - any additives (e.g. interacting proteins) are added from reaction start.
# POST-transcriptional - some additives may be added LATER, after reaction was already going for some time.
# The below flag allows to create additional experiments - to be launched by user later, after (POST) something was added.
# This flag creates additional blank experiments in expnos 322..
# This flag is currently set to work only with IVTNMR (setup_case_list #0 and #1). 
flag_post_series_addon = 0

### Variables of the main script
#===============================

### Setup cases specification
if setup_case == 'HN_31P_iminos_1H_tocsy_HN_31P_iminos_1H':
	referenceExpts = [15, 14, 16, 12, 13, 15, 14, 16, 12] # should match the number in expSet!
	expSet =         [22, 23, 24, 25, 26, 27, 28, 29, 30] # 15N, 31P, ZGSOFAST, zgwg, TOCSY, 15N, 31P, ZGSOFAST, zgwg
	#expTimes =       [12, 6.5, 7, 3.7, 10.5, 12, 6.5, 7, 3.7] # 10.5 <> 61' TOCSY, if want to increase time-spacing
	expTimes =       [13, 6.5, 7, 3.7, 10.5, 13, 6.5, 7, 3.7] # 10.5 <> 61' TOCSY, if want to increase time-spacing
	if flag_post_series_addon:
		expSet = [x+300 for x in expSet]
	
elif setup_case == '31P_iminos_1H_tocsy_31P_iminos_1H':
	referenceExpts = [14, 16, 12, 13, 14, 16, 12]
	expSet =         [22, 23, 24, 25, 26, 27, 28]
	expTimes =       [6.5, 7, 3.7, 10.5, 6.5, 7, 3.7]
	if flag_post_series_addon:
		expSet = [x+300 for x in expSet]

elif setup_case == 'mapper_many_hsqcs':
	referenceExpts = [15, 15, 14, 15, 15, 12, 15, 15, 16] # should match the number in expSet!
	expSet =         [22, 23, 24, 25, 26, 27, 28, 29, 30] # 15N, 15N, 31P, 15N, 15N, zgwg, 15N, 15N, ZGSOFAST
	expTimes =       [12, 12, 6.5, 12, 12, 3.7, 12, 12, 7]

elif setup_case == 'mapper_many_hsqcs_and_tocsy':
	referenceExpts = [15, 15, 14, 15, 15, 12, 15, 15, 16, 15, 15, 13]
	expSet =         [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]
	expTimes =       [12, 12, 6.5, 12, 12, 3.7, 12, 12, 7, 12, 12, 10.5]
		
elif setup_case == 'SNIC_no_15Nprot':
	referenceExpts = [12, 14, 16]
	expSet =         [22, 23, 24]
	expTimes =       [3, 7, 11]
	
elif setup_case == 'SNIC_with_15Nprot':
	referenceExpts = [12, 14, 15, 16]
	expSet =         [22, 23, 24, 25]
	expTimes =       [3, 7, 13, 11]
	
elif setup_case == '5Pz':
	referenceExpts = [15, 14, 12, 13, 15, 16] # should match the number in expSet!
	expSet =         [22, 23, 24, 25, 26, 27] # 15N, 31P, zgwg, TOCSY, 15N, ZGSOFAST
	expTimes =       [32, 6.5, 3.7, 5, 32, 10]

elif setup_case == 'HN_31P_every_second_1D1H':
	referenceExpts = [15, 14, 12, 15, 14] # should match the number in expSet!
	expSet =         [22, 23, 24, 25, 26] # HN, 31P, 1H, HN, 31P
	expTimes =       [900, 6.5, 3.7, 900, 6.5]

elif setup_case == 'HN_31P_every_second_1D1H_or_sofast':
	referenceExpts = [15, 14, 12, 15, 14, 16] # should match the number in expSet!
	expSet =         [22, 23, 24, 25, 26, 27] # HN, 31P, 1H, HN, 31P, SF
	expTimes =       [70, 6.5, 3.7, 70, 6.5, 10]

	#elif setup_case == '':
#	referenceExpts = []
#	expSet =         []
#	expTimes =       []

copyReferenceExptsAndCleanTitles = 1 # 1=yes; 0=no

seriesStart = expSet[-1]+1 # final expt+1. or can set manually: 33 # w/o hsqcs. first experiment number in created series
# When adding experiments after initial creation - set above copyReferenceExptsAndCleanTitles = 0
#seriesStart = 148 # final expt+. or can set manually: 33 # w/o hsqcs. first experiment number in created series

if seriesStart > 50:
	if CONFIRM("Warning",
		'<html>SeriesStart = <font color="red">'+str(seriesStart)+'</font>'+
		'<br><br> Continue?! (press Cancel to stop the script)</html>') == 0:
		EXIT()

#===============================
### Actual script
#===============================
NAME, expno, procno, CURDIR = CURDATA()


## Save a copy of itself to dataset dir
#======================================
py_filename = os.path.basename(sys.argv[0])
py_bak_path = os.path.join(CURDIR,NAME,(py_filename+'.copy'))
print ('Copying %s to %s ...' % (py_filename,py_bak_path))
shutil.copy(sys.argv[0], py_bak_path)


#======================================
print 'Setup case=', setup_case

XCMD('rep 1', wait = WAIT_TILL_DONE) # to avoid errors when trying to start copying from non-first procno
curdat = CURDATA() # get current dataset
datasetFolder = curdat[3]+"/"+curdat[0]+"/"

## Get confirmation - check that target DATASET is correct
if len(sys.argv) > 1 and sys.argv[1] == 'noconfirm':
	print("noconfirm option specified - running w/o user confirmation")
else:
	if CONFIRM("Warning",
		'<html>The script will execute in: <br><br>'+
		'<font color="red">'+curdat[0]+'</font>'+
		'<br><br> Sure this is your TARGET dataset (data m/b be overwritten)! <br> Continue?!</html>') == 0:
		EXIT()

## If reference expts exist - ask to use the existing ones. Or abort copying.
if copyReferenceExptsAndCleanTitles & os.path.isdir(os.path.join(datasetFolder,str(expSet[0]))):
	if CONFIRM("Warning",
		'<html>Starting experiments '+str(expSet[0])+'-'+str(expSet[-1])+' seem to already exist. <br> Press OK to use the existing ones. Or CANCEL to abort.</html>'):
		copyReferenceExptsAndCleanTitles = 0
	else:
		EXIT()


# Create first time point serie with blank titles
if copyReferenceExptsAndCleanTitles and (len(referenceExpts) != len(expSet)):
	ERRMSG(message = "Number of reference expts and expts in expSet must be the same!", title=None, details=None, modal=0) # modal=1 if want the program to pause here until Close button is pressed.
	EXIT()

if copyReferenceExptsAndCleanTitles:
	for exp in range(0,len(referenceExpts)): # Loop to copy every experiment in the block
		XCMD('Re '+str(referenceExpts[exp]), wait = WAIT_TILL_DONE) # read
		XCMD('wraparam '+str(expSet[exp])+' y', wait = WAIT_TILL_DONE) # copy and OVERWRITE! (only the startins expts are overwritten)
		#print('Copied expno '+str(referenceExpts[exp])+' to expno '+str(expSet[exp]))

		# Blank the title
		titleFilePath = datasetFolder+str(expSet[exp])+'/pdata/1/title'
		f = open(titleFilePath, 'w')
		f.write('')
		f.close()

# Create actual time-series
blockTime = sum(expTimes)
timesToRepeat = int( round(timeToRun/blockTime) )
# Diagnostics output to TopSpin shell:
print("Blocksize " + str(blockTime) + " min")
print("Repeat " + str(timesToRepeat) + " times")
#MSG('Repeat '+str(timesToRepeat)+' times')

for exp in range(0,len(expSet)): # Loop to copy every experiment in the block
	targetExp = seriesStart+exp
	XCMD('Re '+str(expSet[exp]), wait = WAIT_TILL_DONE) # read

	for i in range(0,timesToRepeat): # Loop for the required number of experiment repetitions
			XCMD('wraparam '+str(targetExp+len(expSet)*i), wait = WAIT_TILL_DONE) # copy
#			WR(dataset = None, override = "y") # Can try this instead of XCMD('wraparam') - maybe will be more efficient			


# Make the first experiment in post-series - a zg-sofast (imino-detection)
# > To see if there are moderately-quick changes in RNA at this point
if flag_post_series_addon:
	target_exp = 321
	XCMD('Re '+str(16), wait = WAIT_TILL_DONE) # read
	XCMD('wraparam '+str(target_exp)+' y', wait = WAIT_TILL_DONE) # copy and OVERWRITE! (only the startins expts are overwritten)
	# Blank the title
	titleFilePath = datasetFolder+str(target_exp)+'/pdata/1/title'
	f = open(titleFilePath, 'w')
	f.write('')
	f.close()


# Fore debugging use:
#MSG(VARNAME) # displays in a pop-up window
#print("VARNAME = " + str(VAR))


#def main():		

# This is the standard boilerplate that calls the main() function.
# This is not needed if running from TopSpin
#if __name__ == '__main__':
#		main()
