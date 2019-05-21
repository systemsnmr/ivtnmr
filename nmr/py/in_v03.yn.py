# Script for automatic setup of IVTNMR experiments.
## Requires:
#===============
# - An existing "IVTNMR_template" (experiments_template_dir) dataset directory with template experiments and notes_template.txt file.
# - Experiments need to be already optimized for specific spectrometer - the script only does basic things: atmm, topshim, P90H pulse calibrations.
# User expected to
# - do temperature calibration.
# - set correct P90N in the IVTNMR_template (script doesnt do 15N calibration)
# - set correct P90P (31P pulse) (for example after doing paropt calibration/check inside this script).
#   (usually not necessary, cuz 31P pulse is rather stable).


## Comments / KNOWN ERRORS / BUGS
#===============
# - O1 (1H) is currently fixed to 2818 (empirically) - but can turn on the calibration routine - see H_Hz_calibrated
# - Better do not switch around in TopSpin when the script is running.
#   (script strives to be fool-proof, but not guarantee there :-)
# - Script sometimes is waiting until smth happened with notes.txt.
#   If TS shell is hanging on the "=== Log - creating new IVTNMR setup and basic calibrations"
#   Close GEDIT with notes.txt, then the script should proceed.
# - Helper script autolog.py - shows automation log (stored in DATASETFOLDER/automation.og


# TODO:
#===============
# - Figure out why need to close the 'notes.txt' file before the script can proceed?
# - include XCMD('.ret', WAIT_TILL_DONE) everywhere before RE(expt...) places - so that if user moves around, and is by accident in an .md or .ph mode - the script is less likely to get stuck.
# - Option to turn on-off decoupling in zg-iminos -- when there are no Trp iminos in the protein
# - Figure out how to make programmatically ATMM/ATMA just to save wobblecurve - w/o adjusting the tune-match.
# - Move all hard-coded parts into the header. E.g. create expno_2DHN and etc in the header - to not have them hard-coded in the main().
# - Add an option to NOT DO TUNE-MATCH
#   (e.g. if already done manually! could save some time for 1H)
# DONE:
# + Add an option to START experiments after automatic setup.
# + START_TIME0_EXPTS_AFTER_CALIBRATION - make just a flag 0/1 and define experiments to edit below. (e.g. time0_experiments = source_experiments[3:]


## MAYBE IN FUTURE
#===============
# - a variable/argument to specify experiment list (not to have it hard-coded)?
# - check if any of target expts exist (i.e. list dirs in IVTNMR template from 2 to 20) - compare to current
# - add a setup file - like slist - where store: H_Hz (o1), 15N pulse
# - check pulprog name - when setting calibrated P90H, O1 (instead of just defined expno 14)


## Versions:
#===============
# in_v03 - 2019-05-08 @Yar.Nikolaev - cleaned-up for common use
# see earlier versions in in_v02

from __future__ import division # to avoid errors in float division
from __future__ import with_statement # for writing files
from TopCmds import *
from subprocess import call
import os.path
from os import remove
import math
import datetime
import sys # to get own filename
import shutil # to copy files

## SETTINGS
#=========================
#### User to check:
INCLUDE_2DHN = 1 # activate this if not measuring 2DHN
INCLUDE_PAROPT_31P = 1 # If not skipping paropt - need to manually check the P90P value and enter it into the 1D31P experiment!

#### Rarely changed -- THIS IS MORE FOR DEVELOPERS (not users!)
AUTOMATICALLY_ACQUIRE_TIME0_EXPTS_AFTER_CALIBRATION = 1 # Comment out or set to empty if don't want to auto-acquire these.
USE_EXACT_TUNE_MATCH = 0 # Using exact option takes quite a bit longer. Sometimes easier to do by hand instead and not use 'exact'.
SKIP_SLOW_PARTS_FOR_TESTING = 0 # skips paropt, locktuneshim
USE_OWNLOG_FUNC = 1 # A workaround - own logging function, since TopSpin4 has bugs with Python Logging (duplicates logging objects). Maybe will be solved in future.

O1_1H_SETTING = 'fixed_for_all' # Three options: 'keep_as_in_expno' 'fixed_for_all' 'calibrate_automatically'
O1_1H_FIXED_VALUE = 2818 # used only when 'fixed_for_all' setting above

## Variables
#--------
refDB = -9.0 # six2=-9; # 1H channel power level used during Pulsecal. It seemed pulsecal was more accurate at higher power levels - so the refDB here m/b higher than "standard" power level for 1H.
solvent = "H2O+D2O"
experiments_template_dir = 'IVTNMR_template'
# For definition of source experiments to measure - see main() function


#############################################################################
#############################################################################
## Main script
#############################################################################
#############################################################################

NAME, expno, procno, CURDIR = CURDATA()

### Config of logging
log_file_path = CURDIR+'/'+NAME+'/'+'automation.log'

if USE_OWNLOG_FUNC:
	class ManualLog():
		
		def info(self, message):
			f=open(log_file_path, "a")
			log_string = ("%s %s\n" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),message))
			print message
			f.write(log_string)
			f.close()
		
		def debug(self, message):
			self.info(message)

	log = ManualLog()

else:
	import logging as log
	log.basicConfig(
		filename=log_file_path,
		filemode='a', # overwrite file contents instead of appending
		level=log.DEBUG, # minimal level of severity to keep in log: DEBUG INFO WARNING ERROR CRITICAL
		format='%(asctime)s %(levelname)-10s %(message)s',
		#datefmt='%Y%m%d_%H%M%S' # format w/o dashes, but cant easily add millisec
		)

	## Add log display to the console
	log.getLogger().addHandler(log.StreamHandler())

#####################################

### Use this header to debug small things
#if 1: # debugging part
	#src = os.path.join(CURDIR, experiments_template_dir, 'notes_template.txt')
	#dest = os.path.join(CURDIR,NAME,'notes2.txt')
	## If just copying
	#shutil.copyfile(src, dest)
	## Adding title
	#with file(src, 'r') as original: data = original.read()
	#with file(dest, 'w') as modified: modified.write(("\n== %s\n" % NAME) + data)
	#call(["gedit", dest])	# open

	#for e in [12,13,15,16]:
	#	log.debug("expno = %i" % e)

	#d1 = float( GETPAR('D 1') )
	#PUTPAR('D 1', str(3))
	#d1 = float( GETPAR('D 1') )
	#log.debug("D1 = %s" % d1)
	
	#l_a = 0
	#log.debug("lock_attempt = %s" % str(l_a))
	#l_a += 1
	#log.debug("lock_attempt = %s" % str(l_a))
	
	#XCPR('xau paropt.yn "P 1" 80 10 3', WAIT_TILL_DONE)

#EXIT()

#####################################
def copyTemplate(sExpts,tExpts,sSet,tSet):
	log.info('= Copying expts..')
	""" Copies (sources expts) to (target expts) from (template set) to (target set).
	"""
	# TODO - check that target set exists!
	for i in range(0,len(sExpts)): 
		#print(sourceExpts[i])
		source = [sSet, str( sExpts[i] ), "1", CURDIR]
		target = [tSet, str( tExpts[i] ), "1", CURDIR]

		RE(source, show="n") # read the dataset w/o showing it. Syntax/defaults: RE(dataset = None, show = "y")

		# KEEP IN MIND - WR copies also the data!! (not only the parameters)
		# Unlike other functions, WR() gives error if override parameter in not just <"n"> but <override="n">
		# TODO - check if can add WRAPARAM in TopCmds?
		#	WR(target, "n") # write/copy the current dataset. Ask to override if exists. Syntax/defaults: WR(Dataset = None, override = "y") 
		WR(target) # overwrites by default
		# Example of doing ternary if in one line:
		# 'Yes' if fruit == 'Apple' else 'No'
		log.debug('Copied experiment #'+str(i+1)+"\n"+'  '.join(target))

	log.info('Finished copying expts\n')


def readExpt(n,e,p,dir):
	fullpath = [n, str(e), str(p), dir]
	RE(fullpath, show="y") # Syntax/defaults: RE(dataset = None, show = "y")
	log.debug('RE experiment: \n'+'  '.join(fullpath)+'\n')
	
	
def calibrateP1():
	pn,pe,pp,pd = CURDATA()
	log.debug('Starting P1 calibration in ...')
	log.debug(CURDATA())
	
	exptDB = float( GETPAR('PLdB 1') )
	log.debug("exptDB = %.2f" % exptDB)
	exptW = 10**(exptDB/-10.0) # refDB = -10*math.log10(refW)
	log.debug("exptW = %.4f" % exptW)
	# refDB = -9.0 # six2=-9; six=X; five=X
	refW = 10**(refDB/-10.0) # refDB = -10*math.log10(refW)
	PUTPAR('D 1', str(3)) # to be on safe side when using high power for calibration
	PUTPAR('PLW 1', str(refW))
	PUTPAR('P 1', str(5)) # set p1 to low-us, otherwise pulsecal may fail!
	log.debug('refW for pulsecal = %s' % GETPAR("PLW 1"))
	log.debug('Starting P1 for pulsecal = %s' % GETPAR("P 1"))
	XCPR("xau pulsecal quiet same", WAIT_TILL_DONE)
	readExpt(pn,pe,pp,pd)
	P90H = float( GETPAR("P 1") )
	log.debug('P1 @ %.4f W = %.2f' % (refW, P90H))
	P90H = math.sqrt(refW / exptW) * P90H
	log.debug('P1 @ %.4f W = %.2f' % (exptW, P90H))
	PUTPAR('PLW 1', str(exptW))
	PUTPAR('P 1', str(P90H))
	
	log.info('Finished calibration, P1 = %.2f' % P90H)
	return P90H

def calibrateO1():
	log.debug('Starting O1 calibration...')
#	XCPR("xau o1calib", WAIT_TILL_DONE)
	XCPR("o1calib_silent.yn.au", WAIT_TILL_DONE) # find o1calib_silent.yn.au on the six2
	H_Hz = GETPAR("O1")
	log.info('Finished o1 calibration, H_Hz = %s' % str(H_Hz))
	return H_Hz


def lockTuneShim(expno, solvent):
	target = [NAME,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around
	RE(target, show="y") # in case user has switched elsewhere
	
	# TODO: add options for EXACT and tunebxyz

	## LOCKING
	lock_string = "lock "+solvent
	lock_result = -1
	
	lock_attempt = 1
	while lock_result == -1 and lock_attempt < 10:
		log.info('Locking '+lock_string+' ...')
		result = XCPR(lock_string, WAIT_TILL_DONE)
		lock_result = result.getResult()
		log.debug('lock_attempt = %s, lock_result (-1 if fails) = %s' % (str(lock_attempt), str(lock_result)))
		lock_attempt += 1
		
	if lock_result == -1:
		log.error('Could not find lock. Aborting automation')
		# in case of multi-sample setup - should go to the next sample!
		EXIT()

	## ATMA
	RE(target, show="y") # in case user has switched elsewhere
	# was:
	#log.info('ATMA...')
	#result = XCPR("atma exact", WAIT_TILL_DONE) # FIXBACK - now done below

	log.info('= Tune/match 1H...')
	# FIXBACK - comment out next lines to speed up setup if 1H is already tune
	#tunematch(expno,'f1',100) # tune 31P # store f1@100, f2@200, ..
	if USE_EXACT_TUNE_MATCH:
		tunematch2(expno,'f1',procno_to_store=100) # tune 31P # store f1@100, f2@200, ..
	else:
		tunematch2(expno,'f1',unexact_store_wobb=1)

	log.debug('autophase')
	XCPR("autophase", WAIT_TILL_DONE)
	log.debug('autogain')
	XCPR("autogain", WAIT_TILL_DONE)

	log.info('Shimming...')
	RE(target, show="y") # in case user has switched elsewhere
#	XCPR("topshim ls tunebxyz", WAIT_TILL_DONE)
#	XCPR("topshim ls tuneb", WAIT_TILL_DONE)
	XCPR("topshim ls", WAIT_TILL_DONE)
	log.info('Shimming 2nd time...')
	RE(target, show="y") # in case user has switched elsewhere
	XCPR("topshim ls", WAIT_TILL_DONE)
	
	log.debug('autogain')
	XCPR("autogain", WAIT_TILL_DONE)
### If want to use loopadj - make a SILENT version of it!
### But better d n use it - It optimises for best long-term stability, 
# but not for best lineshape, resolution or homogeneity.
#	log.debug('loopadj')
#	XCPR("loopadj", WAIT_TILL_DONE)
	log.info('lockTuneShim() finished')


# TODO: Earlier used to have a separate tunematch() function.. But since its not used can rename tunematch2 > tunematch.
def tunematch2(expno, channel, procno_to_store=0, unexact_store_wobb=0):
	### Can define multiple OPTIONAL options via corresponding variable their names!
	target = [NAME,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around
	RE(target, show="y") # in case user has switched elsewhere
	
	# Launches atma without EXACT -- 
	# so if reasonably close to the optimum, will not do anything, just store the wobble curve
	if unexact_store_wobb:
		if procno_to_store == 0:
			target_procno = get_next_procno(expno,channel)
		else:
			target_procno = procno_to_store
		atma_command = 'atma %s storeWobb %s' % (channel, str(target_procno))
		
	elif procno_to_store == 0:
		# Do not store result, if procno_to_store is not specified
		atma_command = 'atma %s exact' % (channel)
		
	else:
		atma_command = 'atma %s exact storeWobb %s' % (channel, str(procno_to_store))
		
	log.info(atma_command)
	result = XCPR(atma_command, WAIT_TILL_DONE)
	log.info('tunematch() finished')


### Functions to StoreWobb into first unused procno
def get_immediate_subdirectories(a_dir):
    #print os.listdir(a_dir)
    return [name for name in os.listdir(a_dir)
        if os.path.isdir(os.path.join(a_dir, name))]
    #print name
    
def get_last_storeWobb_procno(procno_dirs,channel):
	#print procno_dirs, channel
	channel_increment = int(channel[-1:])*100 # 100, 200, .. - for f1/f2/..
	#print channel_increment
	procno_dirs = map(int, procno_dirs)
	if channel_increment in procno_dirs: # Check that at least one procno for this channel exists
		last_procno = max([x for x in procno_dirs if channel_increment+100 > x > channel_increment-1])
	else:
		last_procno = channel_increment-1
	return last_procno
	#print procno_dirs, channel_increment, last_procno
	
def get_next_procno(expno, channel):
    pdata_path = os.path.join(CURDIR,NAME,str(expno),'pdata')
    procno_dirs = get_immediate_subdirectories(pdata_path)
    #print pdata_path, procno_dirs
    return int(get_last_storeWobb_procno(procno_dirs, channel))+1 # adds 1 - for the unused

def store_wobb(expno, channel):
	### THIS IS NOT USED ANYMORE!!??
    target = [NAME,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around
    RE(target, show="y") # in case user has switched elsewhere
    
    next_procno = get_next_procno(expno,channel)
    EXIT()
    # !!! XCPR can do ATMA, but not ATMM !!!
    # But there m/b a problem if WAIT_TILL_DONE is not followed upon XCMD
    #XCMD('atmm %s' % channel, WAIT_TILL_DONE)
#    XCPR('atma', WAIT_TILL_DONE)
    #XCMD('wbwr %s' % next_procno, WAIT_TILL_DONE) # save
    #SLEEP(3)
    #XCMD('stop', WAIT_TILL_DONE)
																																				
def recO1andH2Oshim(expno,P90H,H_Hz):
	"""
	Records H2O lineshape - after tap and 360 pulses
	(for shim and o1 checks respectively).
	"""

	target = [NAME,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around

	RE(target, show="y") # in case user has switched elsewhere
		
	PUTPAR('O1', str(H_Hz))
	PUTPAR('P 1', '0.02')
	log.debug('--- Set P1 0.02 on:')
	log.debug(CURDATA())
	try:
		os.remove(os.path.join(CURDIR,NAME,str(expno),'fid'))
	except OSError:
		pass

	# TODO JUST USE ZG(), FT(), PK() here?!
	XCPR('zg', WAIT_TILL_DONE)

	RE(target, show="y") # REREAD - in case user has switched elsewhere
	XCPR('fp', WAIT_TILL_DONE) # combined 'qu zgfp' doesn't work!
	XCPR('apk0', WAIT_TILL_DONE)
	WR([NAME,str(expno),"2",CURDIR]) # overwrites procno 2 (complained if was using 
	
	XCPR('rep 1', WAIT_TILL_DONE) # save for quality checks of shimming
	log.info('Recorded and phased water with tip pulse. Starting 360..')

	P360H = float(P90H)*4
	PUTPAR('P 1', str(P360H)) 
	log.debug('--- Set P1 (360deg)='+str(P360H)+' on dset:')
	log.debug(CURDATA())
	log.debug('Deleted the fid (otherwise may ask user for confirmation):')
	log.debug(os.path.join(CURDIR,NAME,str(expno),'fid'))
	try:
		os.remove(os.path.join(CURDIR,NAME,str(expno),'fid'))
	except OSError:
		pass

	XCPR('zg', WAIT_TILL_DONE)
	RE(target, show="y") # REREAD - in case user has switched elsewhere
	XCPR('fp', WAIT_TILL_DONE)
	log.info('recO1andH2Oshim() finished')
	
	
def recExp(expno,P90H,H_Hz):

	# TODO:
	# RPAR ...
	# If TITLE - just set title and return
	# Else:
  #  - set o1, p1
	#  - run additional commands (if any)
	#  - run zg(), etc

	log.debug('Setting up expno = '+str(expno))
	target = [NAME,str(expno),"1",CURDIR]
	RE(target, show="y") # read the dataset w/o showing it. Syntax/defaults: RE(dataset = None, show = "y")
	PUTPAR('P 1', str(P90H))
	log.debug('--- Set P1 (90deg)='+str(P90H)+' on dset:')
	log.debug(CURDATA())
	PUTPAR('O1', str(H_Hz))
	log.info('Recording experiment # '+str(expno))
	XCPR('zg', WAIT_TILL_DONE)
	
	RE(target, show="y") # read dataset, just in case if users are switching around during run
	log.info('Finished recExp() ')

#def nextSample():
#	#### based on nc.yn.py
#	log.info('Ejecting sample ..')
#	XCPR("ej", WAIT_TILL_DONE)
#	SLEEP(25)   # @M.O.Ebert used 20s.
#	log.info('Injecting new sample ..')
#	XCPR("ij", WAIT_TILL_DONE)
#	log.info('Temperature equilibration (30s) ..')
#	SLEEP(30)
#	log.info('Finished nextSample()')

def copy_and_open_notes_template(name):
	src = os.path.join(CURDIR, experiments_template_dir, 'notes_template.txt')
	dest = os.path.join(CURDIR,name,'notes.txt')
	#shutil.copyfile(src, dest) # If just copying
	## Combining notes_template with dataset title
	with file(src, 'r') as original: data = original.read()
	with file(dest, 'w') as modified: modified.write(("\n## %s\n====================================" % name) + data)
	call(["gedit", dest])	# open


def set_title(dataset,text):
	outpath = dataset[3]+'/'+dataset[0]+'/'+dataset[1]+'/pdata/'+dataset[2]
	f = open(outpath+'/title','w')
	f.write(text)
	f.close()


def get_timetag():
	timevar = time.localtime()
	timetag = '%d%02d%02d_%02d%02d'%(timevar[0],timevar[1],timevar[2],timevar[3],timevar[4])
	return timetag

def copy_self_to_dataset_dir():
	py_filename = os.path.basename(sys.argv[0])
	py_bak_path = os.path.join(CURDIR,NAME,(py_filename+'.copy'))
	log.info('Copying %s to %s ...' % (py_filename,py_bak_path))
	shutil.copy(sys.argv[0], py_bak_path)

#def read_sample_list():
#	dirpath = CURDIR+'/'+NAME
#	sys.path.append(dirpath)
#	#print('dirpath='+dirpath)
#	try:
#		from sample_list import solvent, templateSet, sourceExpts, subtract_to_set_expt_numbers, expset
#	except:
#		MSG("Sample list (sample_list.py) not found or has errors")
#		EXIT()
#	log.debug('Solvent = '+solvent)
#	log.debug('TemplateSet = '+templateSet)
#	log.debug('sourceExpts = '+str(sourceExpts))
#	log.debug('expset = '+str(expset))
	
#	return solvent, templateSet, sourceExpts, subtract_to_set_expt_numbers, expset

####################################
# Actual run commands
def main():
	targetSet = NAME
	
	log.info('=== Log - creating new IVTNMR setup and basic calibrations')
	
	copy_self_to_dataset_dir()
	
	if len(sys.argv) > 1 and sys.argv[1] == 'noconfirm':
		print("noconfirm option specified - recording w/o asking for user confirmation")
	else:
		if CONFIRM("Warning",
			'<html>The automation script will execute in: <br><br>'+
			'<font color="red">'+targetSet+'</font>'+
			'<br><br> This should be the TARGET dataset (any data will be overwritten)! <br> Continue?!</html>') == 0:
			EXIT()

	copy_and_open_notes_template(NAME)

## This was used to read setup info from separate file	
#	solvent, templateSet, sourceExpts, subtract_to_set_expt_numbers, expset = read_sample_list()
	templateSet = experiments_template_dir
	if INCLUDE_2DHN:
		sourceExpts = [2,3,4,12,13,14,15,16] # sample name, 1H cal / atmm / topshim, 31P cal, real expts..
	else:
		log.info('= Skipping 15N 2DHN experiment')
		sourceExpts = [2,3,4,12,13,14,16] # sample name, 1H cal / atmm / topshim, 31P cal, real expts..
	
	time0_expnos_to_record = sourceExpts[3:] # skips first three - title/calibration 1H/calibration 31P
	
	# for our "standard" parameter namings - check: edau hsqc15N.all
	BF1 = GETPAR("BF1")
	P90H = 5
	H_Hz = 4.7*float(BF1);
	
	### Copy expts
	sourceExptsFinal = sourceExpts # just a remnant from multi-sample automation
	targetExpts = sourceExptsFinal

	copyTemplate(sourceExptsFinal,targetExpts,templateSet,targetSet)
	#title = '*********  %s  *********' % expset[i][1]
	title = '*****  INXX: 303K pXX, co/post-XXuM XX, TTD77(-S,-RE,+PPase,+DSS,MgCl2)  *****'
	titleExpt = sourceExptsFinal[0]
	set_title([NAME,str(titleExpt),"1",CURDIR], title)
	
	### Do all calibrations in the first expt after title
	log.info('= Lock, tune/match 1H, shim...')
	calibrExpt = targetExpts[1]
	readExpt(targetSet,calibrExpt,"1",CURDIR)
	
	if not SKIP_SLOW_PARTS_FOR_TESTING:
		lockTuneShim(calibrExpt, solvent)

	## Tune/match heteronuclei
	## Syntax: tunematch(expno, channel, procno_to_store_result). procno parameter is optional
	log.info('= Tune/match 31P...')
	if USE_EXACT_TUNE_MATCH:
		#tunematch(4,'f1',100) # tune 31P # store f1@100, f2@200, ..
		tunematch2(4,'f1',procno_to_store=100) # tune 31P # store f1@100, f2@200, ..
	else:
		tunematch2(4,'f1',unexact_store_wobb=1)
			
	if INCLUDE_2DHN:
		log.info('= Tune/match 15N...')
		if USE_EXACT_TUNE_MATCH:
			#tunematch(15,'f3',300) # tune 15N # store f1@100, f2@200, ..
			tunematch2(15,'f3',procno_to_store=300) # tune 15N # store f1@100, f2@200, ..
		else:
			tunematch2(15,'f3',unexact_store_wobb=1)

	log.info('= 1H calibrations...')
	## Read again, just in case user was switching around TopSpin
	readExpt(targetSet,calibrExpt,"1",CURDIR)
	
	H_Hz_default_in_expno = GETPAR("O1")
	
	## Run o1 calibration, catching and suppressing errors if such occur.
	## Normally have it fixed - cuz automatic calibration: takes time, generates errors, introduces phase variability in 1D1H.
	if O1_1H_SETTING=='calibrate_automatically':
		H_Hz_calibrated = None
		while H_Hz_calibrated is None:
			try:
				H_Hz_calibrated = calibrateO1()
			except: # don't care whats the error here
				log.debug('Got exception during o1calib')
				pass
		H_Hz = H_Hz_calibrated
		log.info('o1 calibrated = ' + str(H_Hz))
	elif O1_1H_SETTING=='fixed_for_all':
		H_Hz = O1_1H_FIXED_VALUE
		log.info('o1 (fixed from empirical data)='+str(H_Hz))
	elif O1_1H_SETTING=='keep_as_in_expno':
		H_Hz = H_Hz_default_in_expno # if its not set - won't be used
		log.info('o1 (left as default in expno)='+str(H_Hz))

	log.debug('P90H(before optimization) [us])='+str(P90H))	
	dBH = GETPAR("PLdB 1")
		
	# Read again, just in case user was switching around TopSpin
	####XXX####
	readExpt(targetSet,calibrExpt,"1",CURDIR)
	if not SKIP_SLOW_PARTS_FOR_TESTING:
		P90H = calibrateP1()
	log.debug('P90H(returned)='+str(P90H))
	P90H = GETPAR("P 1")
	log.debug('P90H(after function)='+str(P90H))


	if not SKIP_SLOW_PARTS_FOR_TESTING:
		log.info('= Record o1/shim quality checks')
		recO1andH2Oshim(calibrExpt, P90H, H_Hz)
	# TODO : in recO1.. can add +/- 0.2us expts - to really see if p1 is at optimum
		XCPR("autoshim off", WAIT_TILL_DONE) # did n test this inside of this script!
	
	### Check heteronuclei stability before doing PAROPT
	## d n work yet! (XCMD dis-synchronization issue)
	#store_wobb(4,'f1') # 31P # store_wobb(expno, channel) - will save into next unused procno
	log.info('= Check 31P atmm before paropt...')
	tunematch2(4,'f1',unexact_store_wobb=1) # 31P # store_wobb(expno, channel) - will save into next unused procno

	### 31P signal check (and phase determination)
	# Can put to separate function? set P1 and record expt?
	#### TEST THIS!!!
	target = [NAME,str(4),"1",CURDIR] # save so that can re-read if ppl are switching around  
	RE(target, show="y") # REREAD
	PUTPAR('P 1', str(45))
	log.debug('--- Set 31P P90 = 45us on:')
	log.debug(CURDATA())

	XCPR('zg', WAIT_TILL_DONE)
	RE(target, show="y") # REREAD - in case user has switched elsewhere
	XCPR('efp', WAIT_TILL_DONE) # combined 'qu zgfp' doesn't work!
	
	# Region for aNTP (or gNTP?)
	PUTPAR('F1P', str(-2.88))
	PUTPAR('F2P', str(-7.26))
	XCPR('apk0', WAIT_TILL_DONE)
	XCPR('wrp 2', WAIT_TILL_DONE)

	if INCLUDE_PAROPT_31P:		
		log.info('= Paropt for 31P pulse check...')
		readExpt(targetSet,"4","1",CURDIR)
		initial_P90P = GETPAR("P 1")

		if not SKIP_SLOW_PARTS_FOR_TESTING:
			XCPR('xau paropt.yn "P 1" 80 10 3', WAIT_TILL_DONE)
		log.info('paropt finished')

		# Return the default / initial P90P value
		readExpt(targetSet,"4","1",CURDIR)
		PUTPAR('P 1', initial_P90P)

	else:
		log.info('paropt skipped')

				
	### Set calibrations on actual experiments
	expnos_to_calibrate = sourceExpts[3:]
	for e in expnos_to_calibrate:
		# TODO: make it check pulprog name: pulprog = GETPAR("pulprog") ? then can just go over whole expt set
		target = [NAME,str(e),"1",CURDIR]
		RE(target, show="y") # read the dataset w/o showing it. Syntax/defaults: RE(dataset = None, show = "y")
		if e != 14:
			PUTPAR('P 1', str(P90H))
			PUTPAR('O1', str(H_Hz))
			log.debug('--- Set P90H(p1)=%s & O1=%s on expno %i' % (str(P90H), str(H_Hz), e))
		else:
			PUTPAR('P 3', str(P90H))
			PUTPAR('O2', str(H_Hz))
			log.debug('--- Set P90H(p3)=%s & O2=%s on expno %i (1D31P)' % (str(P90H), str(H_Hz), e))

	# Copy 15N P90N pulse from 2DHN to imino experiment (for decoupling)
	if INCLUDE_2DHN:
		expno_2DHN = 15
		expno_iminos = 16
		RE([NAME,str(expno_2DHN),"1",CURDIR], show="y") # read 2DHN
		P90N = GETPAR("P 21")
		log.debug('- Set P90N=%s to iminos expno %i' % (str(P90N), expno_iminos))
		RE([NAME,str(expno_iminos),"1",CURDIR], show="y") # read imino
		PUTPAR('P 5', str(P90N))


	## Check STABILITY of tune-match after calibrations.
	#===================================================================================
	## Just storing data - does not work(yet)! (due to XCMD dis-synchronization issue).
	## For the moment - using a "hack-around" with unexact_store_wobb=1.
	## In this case if only slightly detuned - atma will not adjust anything.
	#store_wobb(4,'f1') # 31P # store_wobb(expno, channel) - will save into next unused procno
	#store_wobb(15,'f3') # 15N # store_wobb(expno, channel) - will save into next unused procno
	log.info('= Check 1H atmm...')
	tunematch2(3,'f1',unexact_store_wobb=1) # 1H
	
	log.info('= Check 31P atmm...')
	tunematch2(4,'f1',unexact_store_wobb=1) # 31P
	
	if INCLUDE_2DHN:
		log.info('= Check 15N atmm...')
		tunematch2(15,'f3',unexact_store_wobb=1) # 15N # tunematch2(expno, channel,u_s_w=1) - will save into next unused procno

	## Record reference experiments
	#===================================================================================
	if AUTOMATICALLY_ACQUIRE_TIME0_EXPTS_AFTER_CALIBRATION:
		for expno in time0_expnos_to_record:
			if not INCLUDE_2DHN and expno==15:
				pass
			else:
				log.info('Recording experiment # '+str(expno)+' ...')
				RE([NAME,str(expno),"1",CURDIR], show="y")
				XCPR('zg', WAIT_TILL_DONE)
		
		## Check STABILITY of tune-match after time0 reference experiments were recorded.
		#===================================================================================
		log.info('= Check 1H atmm...')
		tunematch2(3,'f1',unexact_store_wobb=1) # 1H
	
		log.info('= Check 31P atmm...')
		tunematch2(4,'f1',unexact_store_wobb=1) # 31P
	
		if INCLUDE_2DHN:
			log.info('= Check 15N atmm...')
			tunematch2(15,'f3',unexact_store_wobb=1) # 15N # tunematch2(expno, channel,u_s_w=1) - will save into next unused procno

	log.info('Automation script finished.\n\n')

### Own checks for getResult
### getResult returns DIFFERENT THINGS - depending on what 
### function is called. I.e. watchout if implementing this!
#	while result.getResult() < 0:
#		print('result in while loop = '+str(result.getResult()))
#		SLEEP(5)

#####################
# Fore debugging use constructions like:
#MSG(VARNAME) # displays in a pop-up window
#print("VARNAME = " + str(VAR))
#EXIT() # make a break-point

#####################
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
		main()
