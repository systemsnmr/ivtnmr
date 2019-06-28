## XXX

## TODO:
#===============
# - check that expnos for atma, shim exist (and also the first_expno_after_addition_of_T7)


## Versions:
#===============
# in_start_01 - 2019-05-17 @Yar.Nikolaev - based on in_v03.yn.py

from __future__ import division # to avoid errors in float division
from __future__ import with_statement # for writing files
from TopCmds import *
from subprocess import call
import os.path
import math
import datetime
import time
import re
import sys # to get own filename

## SETTINGS
#=========================
USE_OWNLOG_FUNC = 1 # A workaround - own logging function, since TopSpin4 has bugs with Python Logging (duplicates logging objects). Maybe will be solved in future.
solvent = "H2O+D2O"
UNEXACT_TUNE_MATCH = 1 # If set to 1 - will adjust only if badly detuned!
lock_expno = 3
atma_1H_expno = 3 # where atma will be executed and wobb curves stored
atma_31P_expno = 4 # where atma will be executed and wobb curves stored
shim_expno = 3 #
first_expno_after_addition_of_T7 = 22 # Not used currently anymore, cuz automatic execution of the first experiment is turned off.

#############################################################################
#############################################################################
## Main script
#############################################################################
#############################################################################

NAME, expno, procno, CURDIR = CURDATA()
target_set = NAME # replacing, cuz NAME too ambigous - can be used for something else too.

## Config of logging
#====================
log_file_path = CURDIR+'/'+target_set+'/'+'automation.log'

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
def readExpt(n,e,p,dir):
	fullpath = [n, str(e), str(p), dir]
	RE(fullpath, show="y") # Syntax/defaults: RE(dataset = None, show = "y")
	log.debug('RE experiment: \n'+'  '.join(fullpath)+'\n')

def lock(expno, solvent):
	target = [target_set,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around
	RE(target, show="y") # in case user has switched elsewhere
	
	## LOCKING
	lock_string = "lock "+solvent
	lock_result = -1 # initial state - will lock until this one changes!

	lock_attempt = 1
	while lock_result == -1 and lock_attempt < 5:
		log.info('Locking '+lock_string+' ...')
		result = XCPR(lock_string, WAIT_TILL_DONE)
		lock_result = result.getResult()
		log.debug('lock_attempt = %s, lock_result (-1 if fails) = %s' % (str(lock_attempt), str(lock_result)))
		lock_attempt += 1
		
	if lock_result == -1:
		log.error('Could not find lock Aborting automation')
		# in case of multi-sample setup - should go to the next sample!
		EXIT()						
	else:
		log.debug('autophase')
		XCPR("autophase", WAIT_TILL_DONE)
		log.debug('autogain')
		XCPR("autogain", WAIT_TILL_DONE)
				
	log.info('lock() finished')


def shim(expno):
	target = [target_set,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around
	log.info('Shimming...')
	RE(target, show="y") # in case user has switched elsewhere
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
	log.info('shim() finished')

				
def tunematch2(expno, channel, procno_to_store=0, unexact_store_wobb=0):
	### Can define multiple OPTIONAL options via corresponding variable their names!
	target = [target_set,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around
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
	log.info('tunematch2() finished')


## Functions to StoreWobb into first unused procno
#=================================================
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
    pdata_path = os.path.join(CURDIR,target_set,str(expno),'pdata')
    procno_dirs = get_immediate_subdirectories(pdata_path)
    #print pdata_path, procno_dirs
    return int(get_last_storeWobb_procno(procno_dirs, channel))+1 # adds 1 - for the unused

#=================================================

def get_time0():
    notes_path = os.path.join(CURDIR,target_set,'notes.txt')
    textfile = open(notes_path, 'r')
    filetext = textfile.read()
    textfile.close()
    matches = re.findall("<time0>(.*)</time0>", filetext)

    if len(matches) > 1:
        ERRMSG("Script aborted: more than one time0 tags found.")    
    else:
        print("notes.txt time0 = %s" % matches[0])
        
        time_split = matches[0].split()
        time0 = datetime.datetime.strptime(' '.join(time_split[:2]), '%Y-%m-%d %H:%M:%S')
        #print time0
    return time0


## Main procedure
#=================================================
def main():
	
	log.info('=== Starting experiment series after addition of T7 ...')
	
	if len(sys.argv) > 1 and sys.argv[1] == 'noconfirm':
		print("noconfirm option specified - recording w/o asking for user confirmation")
	else:
		if CONFIRM("Warning",
			'<html>The automation script will execute in: <br><br>'+
			'<font color="red">'+target_set+'</font>'+
			'<br><br> This should be the TARGET dataset (data will be overwritten)! <br> Continue?!</html>') == 0:
			EXIT()
	
	lock(lock_expno, solvent)
	
	if UNEXACT_TUNE_MATCH:
		### 2019-06-11 - removed ,procno_to_store=107 from all below
		tunematch2(atma_31P_expno,'f1',unexact_store_wobb=1) # tune 31P # store f1@100, f2@200, ..
		tunematch2(atma_1H_expno,'f1',unexact_store_wobb=1) # tune 1H # store f1@100, f2@200, ..
	else:
		tunematch2(atma_31P_expno,'f1') # tune 31P # store f1@100, f2@200, ..
		tunematch2(atma_1H_expno,'f1') # tune 1H # store f1@100, f2@200, ..
		
	shim(shim_expno)
	
	# Check if 5 min (standard dead-time) already passed
	time0 = get_time0()
	tdelta = datetime.datetime.now()-time0
	log.info('Finished lock-tune-shim after %s min (%s seconds)' % (str(tdelta.seconds/60),str(tdelta.seconds)))
	if tdelta.seconds >= 5*60:
		log.info('More than 5min after T7 addition - can launch qumulti DIRECTLY!')
	else:
		seconds_till_start = 5*60-tdelta.seconds
		log.info('Still  %s  seconds until 5min (default time to qumulti the expts).' % str(seconds_till_start))
		## Make "if 1" below if want to use auto-start of the first expno.
		if 0:
			log.info('Waiting until 5min dead time after T7 addition ...')
			SLEEP(seconds_till_start)
			log.info('Starting expts ...')
	#print time0, time_now, tdelta

	if os.path.exists(os.path.join(CURDIR, target_set, str(first_expno_after_addition_of_T7))):
		RE([target_set,str(first_expno_after_addition_of_T7),"1",CURDIR], show="y")
		### Enable this one if want auto-start first expno!
		### Disabled for now, because of TopSpin "qumulti over running queue" error.
		#log.info('Recording experiment # '+str(first_expno_after_addition_of_T7)+' ...')
		#XCPR('zg', WAIT_TILL_DONE)
	else:
		err_msg = ('Expno %s does not exist. Forgot to run createExpSeries? Aborting...' % first_expno_after_addition_of_T7)
		log.info(err_msg)
		ERRMSG(err_msg)

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
