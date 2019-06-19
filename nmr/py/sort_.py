## Script for sorting of time-resolved NMR experiments. 
## Writes the time from start of the experiment into the title file.
## To be executed from within a dataset (not EDTE, etc).
## Can also generate cara repository with spectra links.
## Check that all used pulse-program names are listed in the script!! (defines expt type by pulprog name)

#tt2 - could also invoke inproc2 after execution.

# TODO:
# - Redefine python_dir -- use TopCmds or something to be not dependent on the path!

# - Make script to save itself
# - Add flags for automatic archiving before and deletion of 12--500 dirs after?
# - DO NOT SKIP EXPTS 12-19 even if they are empty! (in case if skipped refs)
#   (e.g. as when non-recording 17-18
#-----
# - LOG the results of sorting!
# - Better check user input of arguments
# - Add OPTIONS parsing

## @Yar.Nikolaev. 2012-..

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
# 2019-05-xx
#   # - multiple improvements: reads expt time from audita.txt; takes into account UTC of EPU/blade, ...
#==================================================
try: # Try-catch to allow tests in normal Python (outside of TopSPin)
    from TopCmds import *
    RUNNING_IN_TOPSPIN = 1
except:
    RUNNING_IN_TOPSPIN = 0
    pass

from subprocess import call
#import sys, os, time, shutil, commands, collections
import os
import time
import datetime
import re
from datetime import timedelta

flag_sort = 1
flag_generate_cara_repo = 0
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

def convert_datetime_to_UTC(t):
    ret = datetime.datetime.strptime(t[0:19],'%Y-%m-%d %H:%M:%S')
    if t[20]=='+':
        ret-=timedelta(hours=int(t[21:23]),minutes=int(t[23:]))
    elif t[20]=='-':
        ret+=timedelta(hours=int(t[21:23]),minutes=int(t[23:]))
    return ret


def read_audita_times(expno_dir):
    # print expno_dir
    fname='%s/audita.txt'%expno_dir
    f = open(fname,'r')

    found_start = 0
    found_end = 0

    # Go line by line - trying to find end-time and start-time (usually end is first)
    for line in f:
        # print line
        
        ## TODO: might want to implement from MMayzel - check for "completed at"
        
        # Won't do this if already found end
        if not found_end:
            if "(   1,<" in line:
                # print 'found end'
                found_end = 1
                end=line.split(',')[1]
                end=end[1:-1] # drop flanking brackets
                end=end.split()[:-1] # convert to array, drop last element (timezone) // cuz TS python cannot get pytz/tzinfo
                ## For TopSpin3.2 - need to remove microseconds, old python, does not recognize:
                end=' '.join(end)[:-4] # [:-4] - drops microseconds!
                # end = datetime.datetime.strptime(end, '%Y-%m-%d %H:%M:%S.%f')
                end = datetime.datetime.strptime(end, '%Y-%m-%d %H:%M:%S')
                # print end
            continue
            
        # Won't do this if already found the start
        if not found_start:            
            if "started at" in line:
                # print 'found start'
                found_start = 1
                start=line.split()[2:4]
                ## For TopSpin3.2 - need to remove microseconds, old python, does not recognize:
                start=' '.join(start)[:-4] # [:-4] - drops microseconds!
                # start = datetime.datetime.strptime(start, '%Y-%m-%d %H:%M:%S.%f')
                start = datetime.datetime.strptime(start, '%Y-%m-%d %H:%M:%S')
                # print start
            continue
            
        if "acquisition in progress" in line:
            start=0
            end=0
            break            

    f.close()

    center = start + (end-start)/2
    return start, center, end    

def get_time0():
    notes_path = os.path.join(CURDIR,NAME,'notes.txt')
    textfile = open(notes_path, 'r')
    filetext = textfile.read()
    textfile.close()
    matches = re.findall("<time0>(.*)</time0>", filetext)

    if len(matches) > 1:
        ERRMSG("Script aborted: more than one time0 tags found.")
    else:
        # matches[0] = '2018-12-21 22:30:12 +0200' # DEBUG        
        print("notes.txt time0 = %s" % matches[0])

        time0_notes_txt = matches[0]
        time0_parts = time0_notes_txt.split(' ')
        n_parts = len(time0_parts)
        time0_field_length = len(time0_notes_txt)

        if n_parts==2 and time0_field_length==19:
            msg_out = """
            Time0 defined w/o timezone.
            Assuming time0 in the notes.txt matches the timezone of audita.txt.
            """
            print(msg_out)
            time0 = datetime.datetime.strptime(' '.join(time0_parts[:2]), '%Y-%m-%d %H:%M:%S')

        elif n_parts==3 and time0_field_length==25 and len(time0_parts[2])==5:
            msg_out = """
            Time0 defined with timezone. 
            Assuming >=TopSpin4 and >=NEO console, and
            converting time0 to UTC to match audita format.
            """
            print(msg_out)
            time0 = convert_datetime_to_UTC(time0_notes_txt)
        else:
            error('Unrecognized time0 format (expecting: YYYY-MM-DD HH:MM:SS x0000). Aborting..\n')

    print "Final time0: " + str(time0)
    # exit() # DEBUG
    return time0

def get_time_difference(expno_dir,time0):
    """
    Calculates time center point of the experiment (from audita.txt).
    Returns timedelta between the center of expt and time0.
    This should result in ~identical time as with getNMRtime Matlab.
    """
    # ideas: stackoverflow - how-to-get-file-creation-modification-date-times-in-python; + stuff from MMayzel

    try:
        _, time_center, _ = read_audita_times(expno_dir)
        timeDiff = time_center-time0

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
        
    except:
        print ('Some Error in audita-time (e.g. acquisition still in progress), skipping')
        return str(0)        

def timeToTitle(experiment,time0):
    """
    Writes the time difference from the start to current experiment into the title file.
    Does not overwrite if the file exists.
    """
    	
    experiment = str(experiment)
    # acqusFilePath = datasetFolder+experiment+'/'+'acqus'
    expno_dir = os.path.join(CURDIR,NAME,experiment)	
    titleFilePath = datasetFolder+experiment+'/pdata/1/title'
    
    # timeDiff = getTimeDifference(acqusFilePath,time0)
    timeDiff = get_time_difference(expno_dir,time0)
    
    f = open(titleFilePath, 'w')
    f.write( timeDiff )
    f.close()
    
    print 'expno', experiment, 'time =', timeDiff
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

### Constructs for debugging
# print 'breakpoint reached', start                
# EXIT() if RUNNING_IN_TOPSPIN else exit()                            
    
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
	main()