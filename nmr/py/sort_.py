## Script for sorting of time-resolved NMR experiments. 
## Writes the time from start of the experiment into the title file.
## To be executed from within a dataset (not EDTE, etc).
## Can also generate cara repository with spectra links.

#tt2 - could also invoke inproc2 after execution.

################
#### If modified pulse-program names -- add those to the list!!
################

# TODO:
# - Implement reading time from audita.txt !!! (cannot rely on timestamp of file modification!!)

####################################
### TODO-Plan for integration with inproc.py (see TopSpin python scripts - Task 24)
# - Make a sort_ivtnmr.py function out of tt.py
# - Spin out into MAIN FUNC
# - Take 3 input args: dset name + sorting range
# - DO NOT SKIP EXPTS 12-19 even if they are empty! (in case if skipped refs)
# - Read the expno in target dataset
# - Execute the rest
# 
# - CALL THE SORT function from inproc - and wait until its finishes?!
####################################

####################################
# TODO - do not skip expts from 12-19 (refs) even if they are empty
# (e.g. as when non-recording 17-18
## @Yar.Nikolaev. 2012-..

## Options of the script:
#==================================================
## ..

## TODO
#==================================================
# - LOG the results of sorting!
# + Use default range? (e.g. 12-300) And set optional using an argument?
# - Better check user input of arguments
# - Functionalize properly - using __main__
# - Add OPTIONS parsing

## Versions
#==================================================
# 2016-09-09 and onwards: Towards proper automation
#  - Automatically reading time0
#  - Reading of expt range as arguments
#  - Removed the 'copyTimeZeroRefs' part
#    (now these are considered automatically as part of sorting range)
#  - Included MSG/printout error if pulseprogram is not in the defined list
#  - (2017-04-28): adds default range of target expts (12-500) if not defined by user

# 2015-05-11: Extended the script:
#  (a) Use different pulprog names for same experiment (now checks NAME in LIST)
#  (b) Ability to skip experiments if some numbers are missing
#  (c) Have to manually adjust IMINO START (if using several imino-detection spectra with diff exch times)
#  (d) Fixed the bug with missing DAYs when calculating TIME for the header!
#		  (see tt_test_time.py for debugging)

# 2015-09-04. @YN:
	# (a) Included option for automatic copying of time=0 expts
	# (b) Included skipping if experiment does not contain acquired data
	
# 2018-09-17. @YN:
	# - added automatic generation of an empty cara repository for 2DHN spectra
#==================================================

from subprocess import call
#import sys, os, time, shutil, commands, collections
import os
import time
import datetime
import re
from datetime import timedelta

#MSG('- Mind that ForkLift when COPYING TOP DIR, not just folders! may change time from the magnet!!')

flag_sort = 1
flag_generate_cara_repo = 1
flag_auto_process_data = 0 # will execute inproc2.py in the end

start_time = time.time()
[NAME, expno, procno, CURDIR] = CURDATA()

####### User-defined parameters for sorting experiments #######

## Adding +1 at the end to get full range.
## Script will automatically skip, if some expnos are missing or dont contain acquired data!
if len(sys.argv) != 2+1:
    #ERRMSG("Error: Script expects 3 arguments - expnos for 1.Target; 2.Ligand+Target; 3.Ligand")
    print(sys.argv)
    default_start = 12
    default_end = 500
    print ('expRange not defined, using default range: %i-%i' % (default_start, default_end))
    expRange = range(default_start,default_end+1)
else:
	expRange = range(int(sys.argv[1]),int(sys.argv[2])+1)
print expRange[0], '..', expRange[-1]
#EXIT()
    
## Lists of possible pulprog names for each expt type
zgPulprogs = ["zg-wg001.yn", "zg-wg.eth", "zg-wg3919-dec.eth", "zgesgp001.yn"]
tocsyPulprogs = ["stocsy003.yn", "stocsy-v.yn", "stocsy002.yn"]
hsqcPulprogs = ["hsqc15N.all", "sfhmqc01.yn", "hsqcetgpsi_02.yn"]
phosPulprogs = ["zgig002.yn", "zgig"]
iminoPulprogs = ["zg-sofast006.yn", "WESF_008g.yn", "zg-sofast004.yn", "IBS_1D_HETSOFAST_007e.yn", "IBS_1D_HETSOFAST_008ec.yn"]

## Set the starting expnos to where the series will be copied
zgStart, tocsyStart, hsqcStart, phosStart, iminoStart = 2000, 3000, 4000, 5000, 6000

writeTitles = True # set False if times have already been written &/or you dont want to change the title. True otherwise.

################################################################################
## Settings for generation of the cara repository
################################################################################
firstSpectrumId = 1
peaklist_id = 1
threshold = 20000 # 20000 - for new TopSpin4 data. Old settings: TOCSY ~50-200. HSQC ~7500 or 20000 (IN44).


################################################################################
####### Actual script #######
################################################################################
curdat = CURDATA() # get current dataset
#datasetFolder = '/Volumes/Data/yar/_eth2/data_NMR/spectra/test/'
datasetFolder = curdat[3]+"/"+curdat[0]+"/"

def sortFolders(expRange):
    directory_check(NAME)
    time0 = get_time0()
    
    ## Checking if target directories for sorting already exist.
    if os.path.exists(datasetFolder+"/"+str(zgStart)):
        ERRMSG("Target directories already exist - cannot overwrite.")
        sys.exit("Target directories already exist - cannot overwrite.")
        EXIT()
	
    ## Initialize initial counters for every experiment type
    zgCount = 0
    tocsyCount = 0
    hsqcCount = 0
    phosCount = 0
    iminoCount = 0
	
    for exp in expRange:
        # If the source experiment does not exist - skip it:
        if not os.path.exists(datasetFolder+"/"+str(exp)):
            print 'expno', exp, 'does not exist, skipping'
            continue

		# If experiment does not have acquired data - skip it:
        if not (os.path.exists(datasetFolder+"/"+str(exp)+"/fid") or os.path.exists(datasetFolder+"/"+str(exp)+"/ser")):
			print 'expno', exp, 'does not have FID or SER data, skipping'
			continue

        exp = str(exp)
        expFolder = datasetFolder + exp		
        XCMD('Re '+exp, wait = WAIT_TILL_DONE)
        
        pulprog = GETPAR("PULPROG")

        # If pulseprogram is not defined (just a ?"spacer" experiment) - skip it:
        if (pulprog == ''):
        	print 'expno', exp, 'PULPROG empty, continue'
        	continue
	
        # Write the time to the title (should happen only for actual experiments, not for "spacers" w/o pulprog)
        if writeTitles:
        	timeToTitle(exp,time0)
														
        if (pulprog in zgPulprogs):
        	# Copies, preserving all modification dates
        	call(["cp", "-rpf", expFolder, datasetFolder+str(zgStart+zgCount)])
        	zgCount += 1
	
        elif (pulprog in tocsyPulprogs):
        	# Copies, preserving all modification dates
        	call(["cp", "-rpf", expFolder, datasetFolder+str(tocsyStart+tocsyCount)])
        	tocsyCount += 1
	
        elif (pulprog in hsqcPulprogs):
        	# Copies, preserving all modification dates
        	call(["cp", "-rpf", expFolder, datasetFolder+str(hsqcStart+hsqcCount)])
        	hsqcCount += 1

        elif (pulprog in phosPulprogs):
        	# Copies, preserving all modification dates
        	call(["cp", "-rpf", expFolder, datasetFolder+str(phosStart+phosCount)])
        	phosCount += 1

        elif (pulprog in iminoPulprogs):
        	# Copies, preserving all modification dates
        	call(["cp", "-rpf", expFolder, datasetFolder+str(iminoStart+iminoCount)])
        	iminoCount += 1

        # If pulseprogram is not among known - alert user and skip it:
        else:
        	print 'expno', exp, 'PULPROG is not known, skipping'
        	MSG('Expno ' + str(exp) + ', ' + pulprog + ' - pulprog uknown - skipping')

def directory_check(NAME):
    if CONFIRM("Warning",
        '<html>The automation script will execute in: <br><br>'+
        '<font color="red">'+NAME+'</font>'+
        '<br><br>Some data may be overwritten! <br> Continue?!</html>') == 0:
        EXIT()


def get_time0():
    notes_path = os.path.join(CURDIR,NAME,'notes.txt')
    textfile = open(notes_path, 'r')
    filetext = textfile.read()
    textfile.close()
    matches = re.findall("<time0>(.*)</time0>", filetext)

    if len(matches) > 1:
        ERRMSG("Script aborted: more than one time0 tags found.")
    else:
        print("notes.txt time0 = %s" % matches[0])
        time0_parts = matches[0].split(' ')
        n_parts = len(time0_parts)

        if n_parts == 3:
            print matches[0] + ': time0 has timezone'
            print 'Assuming Op.System time and time0 in notes.txt were in same timezone - no need to change.'
            print 'TODO: include UTC/timezone awareness - if using time from audita.txt in future.'
            time0 = datetime.datetime.strptime(' '.join(time0_parts[:2]), '%Y-%m-%d %H:%M:%S')
            # print timezone_offset

        elif n_parts == 2:
            print matches[0] + ': time0 has no timezone'

            aug_1_2018 = datetime.date(2018, 8, 1)
            april_1_2019 = datetime.date(2019, 4, 1)

            time0 = datetime.datetime.strptime(matches[0], '%Y-%m-%d %H:%M:%S')

            if time0.date() < aug_1_2018:
                print 'Date before Aug 2018. Assuming audita.txt, Op.System time and time0_notes are same timezone.'
                # Leaves time0 intact

            elif time0.date() >= april_1_2019:
                ERRMSG("Aborting .. (After April 1 2019, expecting time0 in notes.txt to have timezone included!)")

            else:
                print 'Date between Aug 2018 and April 2019. UTC difference not fixed. Assuming time0_notes is in UTC (like in EPU/blade), and system time (real time) is one hour later - +0100 (CET).'
                time0 = time0 + timedelta(hours=1)
                print 'time0 after +1-hour correction:' + str(time0)

        else:
            ERRMSG("Something went wrong with parsing time0 stamp, aborting..")
    
    print "Final time0: " + str(time0)
    return time0


def getTimeDifference(filename,time0):
	"""
	Calculates time difference between start of time-series (time0) and end time of current NMR experiment.
	Returns the difference as a _string_ of the form Minutes:Seconds.
	
	Times calculated here provide only a rough timepoint.
	A more accurate way to calculate time-point is to take acquisition start & end times from TopSpin acquisition
	file, and then calculate a middle point from that. This is done in MatLab getNMRdata routine. 
	"""
	# some part taken from http://stackoverflow.com/questions/237079/how-to-get-file-creation-modification-date-times-in-python
	t = os.path.getmtime(filename)
	timeAtExperimentEnd = datetime.datetime.fromtimestamp(t)
    # print timeAtExperimentEnd
	timeDiff = timeAtExperimentEnd-time0
	
	# 2015-05-11: added DAYS in the diffMin calculation:
	diffMin = timeDiff.days*24*60 + (timeDiff.seconds / 60)
	diffSec = timeDiff.seconds - diffMin*60
	#return str(diffMin) + ':' + str(diffSec)
	# Add zero in front, if diffMin is only two-digit number
	if diffMin<100:
		diffMin = "00"+str(diffMin)
	elif diffMin<1000:
		diffMin = "0"+str(diffMin)		
	return str(diffMin)


def timeToTitle(experiment,time0):
	"""
	Writes the time difference from the start to current experiment into the title file.
	Does not overwrite if the file exists.
	"""
	experiment = str(experiment)
	acqusFilePath = datasetFolder+experiment+'/'+'acqus'
	titleFilePath = datasetFolder+experiment+'/pdata/1/title'

	## If title exists - do not overwrite. Disabled for the moment (2016-09-09) - d n make sense anyway!?
    # if os.path.exists(titleFilePath):
	if 0:
		print 'expno', experiment, ': title file exists - will not overwrite.'
	else:
		timeDiff = getTimeDifference(acqusFilePath,time0)
		f = open(titleFilePath, 'w')
		f.write( timeDiff )
		f.close()

	f = open(titleFilePath, 'r')
	print 'expno', experiment, 'time =', f.readline()
	#MSG(timeDiff)


################################################################################
## Generate cara repo
################################################################################

def generateSpectraXML(experimentRange,short_name):
	specId = firstSpectrumId
	out = ""

	for experiment in experimentRange:
		experiment = str(experiment)
		
		titleFilePath = os.path.join(CURDIR,NAME,experiment,'pdata/1/title')
        # print titleFilePath
		
        # titleFilePath = datasetFolder+experiment+'/pdata/1/title'
		f = open(titleFilePath, 'r')
		expName = short_name+'_'+f.readline()
		f.close()
		
        # #TMP-from unsorted
        # if len(experiment) == 2:
        #     expName = short_name+'_'+'0'+experiment
        # else:
        #     expName = short_name+'_'+experiment

        # dataFilePath = datasetFolder+experiment+'/pdata/1/2rr'
		dataFilePath = os.path.join(CURDIR,NAME,experiment,'pdata/1/2rr')
		
#		out = str(specId)+", "+expName+", "+dataFilePath
		
		## 2D TOCSY
#		out += "<spectrum type='2D TOCSY' id='%s' name='%s' path='%s' sample='0'>\n<level pmax='' pnoise='' nmax='' nnoise='' thres='%s' factor='1.000000'/>\n<cal dim='0' off='0' width='0' fold='A'/>\n<cal dim='1' off='0' width='0' fold='A'/>\n</spectrum>\n" % (specId, expName, dataFilePath, threshold)

# Not sure what is this thing doing here!
#		notesPath = "file://%f/%n" % (folderPath, "/notes.txt")

		## 2D 15N-HSQC
		out += "<spectrum type='2D HSQC15N' id='%s' name='%s' path='%s' sample='0'>\n<level pmax='' pnoise='' nmax='' nnoise='' thres='%s' factor='1.000000'/>\n<cal dim='0' off='0' width='0' fold='A'/>\n<cal dim='1' off='0.000000' width='0' fold='A'/>\n</spectrum>\n" % (specId, expName, dataFilePath, threshold)
		
		specId += 1
	#print out
	return out

def generate_peaklist_xml(experimentRange,short_name):
	specId = firstSpectrumId
	datestamp = datetime.datetime.today().strftime('%y%m%d')
	
	out = ""
	out += "<peaklist id='%s' name='%s_%s' home='%s'>\n<dim atom='H1'/>\n<dim atom='N15'/>\n<model id='0' kind='2'>\n<param gain='0.000000' bal='1.000000' width='0.500000' tol='0.000000'/>\n<param gain='0.000000' bal='1.000000' width='0.500000' tol='0.000000'/>\n</model>\n" % (peaklist_id, datestamp, short_name, specId)

	for experiment in experimentRange:
		out += "<spec id='%s'/>\n" % (specId)
		specId += 1
	out += "</peaklist>"
	#print out
	return out
	#print ' '
	print '### Can also add peaks into peaklist - using <peak id=".."> - but need individual template for that'

######################################################

def get_immediate_subdirectories(a_dir):
    """..."""
    # print os.listdir(a_dir)
    return [name for name in os.listdir(a_dir)
        if os.path.isdir(os.path.join(a_dir, name))]
    # print name

def int_filter( someList ):
    """Takes a list (of strings), returns a list of integers.
    Filters out non-integers on the way"""
    for v in someList:
        try:
            x = int(v)
            yield x # Keep these
        except ValueError:
            continue # Skip these

def get_recorded_expnos(curdir,name):
    """Takes system path, returns structure with thousand-based expno numbers.
    In this script we're only using 2DHN - 4xxx"""
    
    all_dirs = get_immediate_subdirectories( os.path.join(curdir,name) )
    # print all_dirs

    dir_numbers = list( int_filter( all_dirs )) # convert to numbers, dropping non-integers
    # print dir_numbers

    expno_sets = []

    for i in range(1,10):
        incr = i*1000
        # print i, incr
        # Selecting not until x999 - cuz there m/b ref expts (e.g. 4998, 4999 for 31P) there
        expno_sets.append([e for e in dir_numbers if incr+900 >= e >= incr])
    # print expno_sets
    return expno_sets

######################################################
######################################################
######################################################

def main():
    #####################
    if flag_sort:        
        print '= Starting the sorting procedure'
        sortFolders(expRange)
        print('Sorting of %s complete' % NAME)
        print('Sorting completed in %.4f sec' % (time.time()-start_time))        
    
    #####################
    if flag_generate_cara_repo:
        print('Generating cara repository ...')    

        short_name = re.findall('IN[A-Za-z0-9]{0,}', NAME)[0]
        # automatically reads expnos - ALL sets - 2xxx, 3xxx, 4xxx
        expno_sets = get_recorded_expnos(CURDIR,NAME) 
        # indexing starts with 0, so 4xxx expts will have index 4-1: need to decrement
        experimentRange = expno_sets[4-1] 

        spectra_xml = generateSpectraXML(experimentRange,short_name)
        peaklist_xml = generate_peaklist_xml(experimentRange,short_name)
        # print spectra_xml
        # print peaklist_xml

        ### Combine into one file    
        # If the file already exists - create a new with current date/time (to not overwrite)
        cara_repo = curdat[3]+"/"+curdat[0]+"/"+short_name+".cara"
        if os.path.isfile(cara_repo):
            datetime_append = datetime.datetime.today().strftime('%y%m%d_%H%M')
            cara_repo = curdat[3]+"/"+curdat[0]+"/"+short_name+"_"+datetime_append+".cara"

        ## Get parts of cara repo from py/user/data_modules
        # assuming:
        # - getcwd returns "TOPSPIN_INSTALLATION_DIR/prog/curdir/username"
        # - python are located in TOPSPIN_INSTALLATION_DIR/exp/stan/nmr/py/user/
        # - generate path to data_modules:

        # removes three last dirs:
        #/opt/topspin3.2/prog/curdir/yar
        # need *list_name statement to use the list
        cwd_elements = os.getcwd().split( os.sep )
        topspin_installation_dir = os.path.join( *cwd_elements[0:-3] )
        python_dir = 'exp/stan/nmr/py/user'
        ## On our servers - /mol/imb/tp3 - issues with filesystem access - cannot read in subdirs.
        ## - putting for now files directly to the py/user
        # modules_dir = os.path.join( os.sep, topspin_installation_dir, python_dir, 'data_modules')
        modules_dir = os.path.join( os.sep, topspin_installation_dir, python_dir)
    
        # write to the file
        # print cara_repo
        f = open(cara_repo, 'w')
        f_top = open( os.path.join( modules_dir, '180917_INx_blank_1top.cara') , 'r')
        f_botoom = open( os.path.join( modules_dir, '180917_INx_blank_2bottom.cara'), 'r')
        f.write( f_top.read() )
        f.write( spectra_xml )
        f.write( peaklist_xml )
        f.write( f_botoom.read() )
        f.close()
        print('...done')  
        
    if flag_auto_process_data:
        EXEC_PYFILE('inproc2.py')  
    
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
	main()