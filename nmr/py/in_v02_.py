# Disclaimer: Limitations of Liability for this code
# This script is provided "as is" - merely as a guideline for automated setup. Its working on our 600 AVANCEIII / NEO consoles. But it was not thoroughly tested on different spectrometers and may contain bugs and incompatibilities with TopSpin versions. Authors assume no responsibility, and shall not be liable to you or to any third party for any direct, indirect, special, consequential, indirect or incidental losses, damages, or expenses, directly or indirectly relating to the use or misuse of the code and pulseprograms provided here.

# Better DO NOT SWITCH AROUND in TOPSPIN while script is running!

# See SKIP_SLOW & UNEXACT_TUNE_MATCH variables below - to tweak how fast the script will run

# ERRORS
# - Check - seems script sometimes is waiting until smth happened with notes.txt - then just close the notes.txt file!
# > If TS shell is hanging on the "=== Log - creating new IVTNMR setup and basic calibrations"
# > Close GEDIT with notes.txt - should proceed.

# This script expects:
# - User to do temperature calibration
# - Correct P90N set in the IVTNMR_template (script d n do 15N calibration)

# TODO:
# - Add copying 15N pulse from hmqc to sofast (decoupling)
# - Add an option to NOT DO TUNE-MATCH
#   (e.g. if already done manually! could save some time for 1H)
#########################

from __future__ import division # to avoid errors in float division
from __future__ import with_statement # for writing files
from TopCmds import *
from subprocess import call
import os.path
from os import remove
import math
import logging as log
#import shutil

MSG('pulsecal doesnt work as expected on the NEO, calibrate p1 manually, and set it to correct value in IVTNMR_template/3 before starting this script').
## - Calibration of P90H is currently commented out - find this section - ####XXX####

### "TOP-LEVEL" variables (technically NOT GLOBAL)

CALIBRATE_ON_PO4 = 0 # KEEP THIS OFF GENERALLY!!!
SKIP_PAROPT_31P = 1
SKIP_SLOW_PARTS_FOR_TESTING = 0 # skips paropt, locktuneshim
SKIP_2DHN = 1 # activate this if not measuring 2DHN
UNEXACT_TUNE_MATCH = 1 # quick hack to NOT do TUNE-MATCH - if consecutive days
#DEBUG = 1 # Not using lately. For code debugging - set to zero for normal runs.

NAME, expno, procno, CURDIR = CURDATA()

### Config of logging

log_file = CURDIR+'/'+NAME+'/'+'automation.log'

log.basicConfig(
	filename=log_file,
	filemode='a', # overwrite file contents instead of appending
	level=log.DEBUG, # minimal level of severity to keep in log: DEBUG INFO WARNING ERROR CRITICAL
	format='%(asctime)s %(levelname)-10s %(message)s',
	#datefmt='%Y%m%d_%H%M%S' # format w/o dashes, but cant easily add millisec
	)

## Add log display to the console
## Duplication in the console is not due to this - something else is causing it (errors?)
log.getLogger().addHandler(log.StreamHandler())

#####################################

### Use this header to debug small things
#if 1: # debugging part
	#src = os.path.join(CURDIR,'IVTNMR_template','notes_template.txt')
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
	log.debug('Starting P1 calibration...')
	
	exptDB = float( GETPAR('PLdB 1') )
	log.debug("exptDB = %.2f" % exptDB)
	exptW = 10**(exptDB/-10.0) # refDB = -10*math.log10(refW)
	log.debug("exptW = %.4f" % exptW)
	refDB = -9.0 # six2=-9; six=X; five=X
	refW = 10**(refDB/-10.0) # refDB = -10*math.log10(refW)
	PUTPAR('D 1', str(3)) # to be on safe side when using high power for calibration
	PUTPAR('PLW 1', str(refW))
	PUTPAR('P 1', str(5)) # set p1 to low-us, otherwise pulsecal may fail!
	log.debug('refW for pulsecal = %s' % GETPAR("PLW 1"))
	log.debug('Starting P1 for pulsecal = %s' % GETPAR("P 1"))
	XCPR("xau pulsecal quiet same", WAIT_TILL_DONE)
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
	log.info('Finished calibration, H_Hz = %s' % str(H_Hz))
	return H_Hz


def lockTuneShim(expno, solvent):
	target = [NAME,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around
	
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
	if not UNEXACT_TUNE_MATCH:
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
	# TODO - set AUTOSHIM ON somehow


def tunematch(expno, channel, procno_to_store=0):
	target = [NAME,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around
	
	## ATMA
	RE(target, show="y") # in case user has switched elsewhere
	# Do not store result, if procno_to_store is not specified
	if procno_to_store == 0:
		atma_command = 'atma %s exact' % (channel)
	else:
		atma_command = 'atma %s exact storeWobb %s' % (channel, str(procno_to_store))		
	log.info(atma_command)
	result = XCPR(atma_command, WAIT_TILL_DONE)

	log.info('tunematch() finished')


#"""
def tunematch2(expno, channel, procno_to_store=0, unexact_store_wobb=0):
	### Can define multiple OPTIONAL options via corresponding variable their names!
	target = [NAME,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around
	RE(target, show="y") # in case user has switched elsewhere
	
	if unexact_store_wobb:
		next_procno = get_next_procno(expno,channel)
		atma_command = 'atma %s storeWobb %s' % (channel, str(next_procno))
	
	elif procno_to_store == 0:
		# Do not store result, if procno_to_store is not specified
		atma_command = 'atma %s exact' % (channel)
		
	else:
		atma_command = 'atma %s exact storeWobb %s' % (channel, str(procno_to_store))
		
	log.info(atma_command)
	result = XCPR(atma_command, WAIT_TILL_DONE)
	log.info('tunematch() finished')
#"""


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
	src = os.path.join(CURDIR,'IVTNMR_template','notes_template.txt')
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
	solvent = "H2O+D2O"
	templateSet = "IVTNMR_template"
	if not SKIP_2DHN:
		sourceExpts = [2,3,4,12,13,14,15,16] # sample name, 1H cal / atmm / topshim, 31P cal, real expts..
	else:
		log.info('= Skipping 15N 2DHN experiment')
		sourceExpts = [2,3,4,12,13,14,16] # sample name, 1H cal / atmm / topshim, 31P cal, real expts..
	
	# for our "standard" parameter namings - check: edau hsqc15N.all
	BF1 = GETPAR("BF1")
	P90H = 5
	H_Hz = 4.7*float(BF1);
	
	### Copy expts
	sourceExptsFinal = sourceExpts # just a remnant from multi-sample automation
	targetExpts = sourceExptsFinal

	copyTemplate(sourceExptsFinal,targetExpts,templateSet,targetSet)
	#title = '*********  %s  *********' % expset[i][1]
	title = '*****  INXX: 303K pXX, co/post-XXuM XX, TTD77(-S,-RE,+PPase,+DSS,24mM MgCl2)  *****'
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
	if not UNEXACT_TUNE_MATCH:
		#tunematch(4,'f1',100) # tune 31P # store f1@100, f2@200, ..
		tunematch2(4,'f1',procno_to_store=100) # tune 31P # store f1@100, f2@200, ..
	else:
		tunematch2(4,'f1',unexact_store_wobb=1)
			
	if not SKIP_2DHN:
		log.info('= Tune/match 15N...')
		if not UNEXACT_TUNE_MATCH:
			#tunematch(15,'f3',300) # tune 15N # store f1@100, f2@200, ..
			tunematch2(15,'f3',procno_to_store=300) # tune 15N # store f1@100, f2@200, ..
		else:
			tunematch2(15,'f3',unexact_store_wobb=1)

	log.info('= 1H calibrations...')
	## Read again, just in case user was switching around TopSpin
	readExpt(targetSet,calibrExpt,"1",CURDIR)
	## Not using this for the moment - o1 is ~stable for given temperature and buffer
##	## Run o1 calibration, catching and suppressing errors if such occur.
##	H_Hz_calibrated = None
##	while H_Hz_calibrated is None:
##		try:
##			H_Hz_calibrated = calibrateO1()
##		except: # don't care whats the error here
##			log.debug('Got exception during o1calib')
##			pass
		
##		H_Hz = H_Hz_calibrated
	H_Hz = 2822
	log.info('o1 (fixed from empirical data)='+str(H_Hz))
	
	log.debug('P90H(before optimization) [us])='+str(P90H))	
	dBH = GETPAR("PLdB 1")
		
	# Read again, just in case user was switching around TopSpin
	readExpt(targetSet,calibrExpt,"1",CURDIR)
	####XXX####if not SKIP_SLOW_PARTS_FOR_TESTING:
		####XXX####P90H = calibrateP1()
	####XXX####log.debug('P90H(returned)='+str(P90H))
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
	if CALIBRATE_ON_PO4:
		# Region for PO4
		PUTPAR('F1P', str(4))
		PUTPAR('F2P', str(0))
	XCPR('apk0', WAIT_TILL_DONE)
	XCPR('wrp 2', WAIT_TILL_DONE)

	if not SKIP_PAROPT_31P:		
		log.info('= Paropt for 31P pulse check...')
		readExpt(targetSet,"4","1",CURDIR)
	
		if not SKIP_SLOW_PARTS_FOR_TESTING:
			XCPR('xau paropt.yn "P 1" 80 10 3', WAIT_TILL_DONE)
		log.info('paropt finished')
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

	### Check heteronuclei stability at the end
	## d n work yet! (XCMD dis-synchronization issue)
	#store_wobb(4,'f1') # 31P # store_wobb(expno, channel) - will save into next unused procno
	#store_wobb(15,'f3') # 15N # store_wobb(expno, channel) - will save into next unused procno
	log.info('= Check 1H atmm...')
	tunematch2(3,'f1',unexact_store_wobb=1) # 1H
	
	log.info('= Check 31P atmm...')
	tunematch2(4,'f1',unexact_store_wobb=1) # 31P
	
	if not SKIP_2DHN:
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
