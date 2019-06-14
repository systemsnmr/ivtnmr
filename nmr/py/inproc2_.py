## Full processing of IVTNMR datasets
## @Yar.Nikolaev. 2017-..

## This script allows to process multiple datasets: if "dsets" variable defined in main()
## See Evn "Topspin Python scripts" for development notes.

#### TODO:
# - Save - TopSpin CORBA error -- fix by defining it LOCAL RUNNING?!

#### Longer-term TODO
# - Make that takes all relevant parameters from first expnos in the series (i.e. not just phc0,etc)
# - Taking DSET names as input user arguments!
# - Take processing settings as an input?

# > Additional "Symptom" for above two: when using function log.info(..) was adding up on each round
# - Logging with the use of separate function does not work!
# - Saving processing.log into a separate file for each dataset (need a separate log-obj for each?)

#### TOPSPIN "BUGS"
# - !!! Sometimes SR is reset to 0. Still could not trace why :()
# - apks phase of 1D31P spectrum is sometimes inverted (signals are negative)

import sys
import os
import logging as log
import time
from subprocess import call # to open log file
import sys # to get own filename
import shutil # to copy files

start_time = time.time()

#####################################
## Settings
#####################################
ref_spectrum_for_2DHN_phase_and_STSI = 4001
ref_spectrum_for_1D1H_phase = 2001
ref_spectrum_for_1D1H_SR = ref_spectrum_for_1D1H_phase

# Can disable some sections - if reprocessing
f_process_1D31P = 1
f_process_iminos = 1
f_process_1D1H = 1
f_process_2DHN = 1

# force_nc_proc = 6 # comment out this variable if don't want to force-fix it in processing.

# Logging level:
log_level=log.DEBUG # Set to log.DEBUG if want to get more details.

#####################################
## Actual script
#####################################

## Constants, dictionaries, etc
#==============================
# IUPAC / from Wishart
DSS_SR_1H15N=0.101329118
DSS_SR_1H31P=0.404808636
DSS_SR_1H13C=0.251449530

bc_mod_dict = {
  "0": "no",
  "1": "single",
  "2": "quad",
  "3": "spol",
  "4": "qpol",
  "5": "sfil",
  "6": "qfil",
}

## Functions
#===========================

# def start_log():
#   """Not used yet.
#   Meant to take target dset dir and keep log of processing there"""
#   #### - uncomment such lines below if want to log into a file
#   name, expno, procno, curdir = CURDATA()
#   log_file = curdir+'/'+name+'/'+'processing.log'
#   log.basicConfig(
#       filename=log_file,
#         # filemode='w', # overwrite file contents instead of appending (default mode is 'a')
#         level=log_level, # minimal level of severity to keep in log: DEBUG INFO WARNING ERROR CRITICAL
#         ## if want to debug something specific - use this level + log.warning('...') in the code
#         # level=log.WARNING,
#       format='%(asctime)s %(levelname)-10s %(message)s',
#       #datefmt='%Y%m%d_%H%M%S' # format w/o dashes, but cant easily add millisec
#       )
#     
#       ## Add log display to the console
#   log.getLogger().addHandler(log.StreamHandler())
# 
#       # log example:
#       # log.debug('Copied experiment #'+str(i+1)+"\n"+'  '.join(target))

def read_expt(n,e,p,dir):
    fullpath = [n, str(e), str(p), dir]
    RE(fullpath, show="y") # Syntax/defaults: RE(dataset = None, show = "y")


def proc_2DHN(expnos, target_procno, sr, proc_params, proc_string):
    """Process series of 2D HN spectra (expnos)"""
    name, expno, procno, curdir = CURDATA() # get current dataset
    # dataset_folder = curdir+"/"+name+"/"    

    sr_1H = float( sr )
    sr_15N = sr_1H*DSS_SR_1H15N
    # print sr_1H, sr_15N

    for e in expnos:
        # Read expt only if it exists (skips empty expts):

        log.info('= Processing expno %s...' % (str(e)))
        # print proc_params['STSI'], proc_params['STSR']
        
        if 1: #os.path.exists(dataset_folder+"/"+str(expno)):
            # XCMD('re '+str(e), wait = WAIT_TILL_DONE) # read
            RE([name, str(e), '1', curdir])
            
            #### SR ####
            # Setting also SF HERE! (because SR is NOT A PARAMETER by itself)
            # PUTPAR("1 SR", str(sr_15N)) # does not work anymore!
            # PUTPAR("SF", GETPAR("BF1")) # resets SR to zero

            sf_channel2 = float(GETPAR('BF2'))+sr_15N/1e6 # SF & BF in MHz, sr in Hz
            PUTPAR("2 SF", str(sf_channel2))            
            sf_channel1 = float(GETPAR('BF1'))+sr_1H/1e6 # SF & BF in MHz, sr in Hz
            PUTPAR("1 SF", str(sf_channel1))
            
            XCMD( ('2 sr %s' % str(sr_1H)), wait = WAIT_TILL_DONE)
            XCMD( ('1 sr %s' % str(sr_15N)), wait = WAIT_TILL_DONE)
            
            if e==expnos[0]: # displays this info only for the first expno            
                log.debug('---- SR ----')
                log.debug('sf1=%s' % str(sf_channel1))
                log.debug('sf2=%s' % str(sf_channel2))
                log.debug('sr1=%s' % str(sr_1H))
                log.debug('sr2=%s' % str(sr_15N))                
                log.debug('------------')

            for param in proc_params:
                PUTPAR( ("2 %s" % param) , str(proc_params[param][0]) )
                PUTPAR( ("1 %s" % param) , str(proc_params[param][1]) )
                
                if e==expnos[0]: # displays this info only for the first expno
                    log.debug("2 %s %s" % (param , str(proc_params[param][0])))
                    log.debug("1 %s %s" % (param , str(proc_params[param][1])))
                
                
                
                
            ## Manually set dimenstion parameters again
            ## otherwise when these are altered, TopSpin sets them properly only on a second round!
            for dim_param in ['SI', 'STSI', 'STSR']:
                PUTPAR( ("2 %s" % dim_param) , str(proc_params[dim_param][0]) )
                PUTPAR( ("1 %s" % dim_param) , str(proc_params[dim_param][1]) )

            ### Do not try to overwrite itself - skip this section if target_procno == 1
            if int(target_procno) != 1:
                # XCMD(('wrp %s' % target_procno), wait = WAIT_TILL_DONE) # write
                WR([name, str(e), str(target_procno), curdir], "y") # overwrite
                # XCMD(('rep %s' % target_procno), wait = WAIT_TILL_DONE) # read            
                RE([name, str(e), str(target_procno), curdir])

            # if 'force_nc_proc' in locals():
            #     # XCMD('xfb n nc_proc -6', wait = WAIT_TILL_DONE) # process and set nc_proc scaling same for all
            #     XCMD(('xfb nc_proc %s n' % force_nc_proc), wait = WAIT_TILL_DONE) # process                
            # else:                
            #     XCMD('xfb n', wait = WAIT_TILL_DONE) # process                
            #     ### If wanted NON-FT-PROCESSED DATA:
            #     # XCMD('xtrf n', wait = WAIT_TILL_DONE) # process

            proc_cmd_array = proc_string.split(';')
            
            for proc_cmd in proc_cmd_array:
                # strip() - removes starting and trailing whitespaces
                XCMD(proc_cmd.strip(), wait = WAIT_TILL_DONE) # process
                ### If wanted NON-FT-PROCESSED DATA:
                # XCMD('xtrf n', wait = WAIT_TILL_DONE) # process
                
            # log.info("Finished proc. NC_proc = %i" % int(GETPAR("NC_proc")))

def proc_1D1H_for_SREF(expno):
    """Process 1D1H - for calibration stuff"""
    XCMD('re '+str(expno), wait = WAIT_TILL_DONE) # read
    name, expno, procno, curdir = CURDATA() # get current dataset
    
    PUTPAR("SI", "128k")
    PUTPAR("BCFW", str(0.2)) # in case if using BC_mod qfil
    PUTPAR("BC_mod", "qfil")
    PUTPAR("WDW", "EM")
    PUTPAR("LB", "0.5")
    # EFP() # using EFP here makes DSS-peak recognition worse!!
    FP()
    APK() # tried apk, apks, apkm - apk seemed most reasonable


def proc_1D(expnos, target_procno, sr, proc_params, proc_string):    
    """Generic function to process 1Ds - 1D31P, 1D1H, iminos"""
    name, expno, procno, curdir = CURDATA() # get current dataset
    # dataset_folder = curdir+"/"+name+"/"    

    for e in expnos:
        # Read expt only if it exists (skips empty expts):

        log.info('= Processing expno %s...' % (str(e)))

        if 1: #os.path.exists(dataset_folder+"/"+str(expno)):
            # XCMD('re '+str(e), wait = WAIT_TILL_DONE) # read
            RE([name, str(e), '1', curdir])

            #### SR ####
            # Setting also SF HERE! (because SR is NOT A PARAMETER by itself)
            # PUTPAR("1 SR", str(sr_15N))
            # PUTPAR("SF", GETPAR("BF1")) # resets SR to zero            
            # Set SF
            sf_channel1 = float(GETPAR('BF1'))+sr/1e6 # SF & BF in MHz, sr in Hz
            PUTPAR("1 SF", str(sf_channel1))
            # Set SR
            XCMD( ('1 sr %s' % str(sr)), wait = WAIT_TILL_DONE)
            
            if e==expnos[0]: # displays this info only for the first expno                
                log.debug('---- SR ----')
                log.debug('sf=%s' % str(sf_channel1))
                log.debug('sr=%s' % str(sr))            
                log.debug('------------')

            for param in proc_params:
                if e==expnos[0]: # displays this info only for the first expno                
                    log.debug("1 %s %s" % (param, str(proc_params[param][0])))
                    
                PUTPAR( ("1 %s" % param) , str(proc_params[param][0]) )

            ### Do not try to overwrite itself - skip this section if target_procno == 1
            if int(target_procno) != 1:
                # XCMD(('wrp %s' % target_procno), wait = WAIT_TILL_DONE) # write
                WR([name, str(e), str(target_procno), curdir], "y") # overwrite
                # XCMD(('rep %s' % target_procno), wait = WAIT_TILL_DONE) # read            
                RE([name, str(e), str(target_procno), curdir])

            proc_cmd_array = proc_string.split(';')
            
            for proc_cmd in proc_cmd_array:
                # strip() - removes starting and trailing whitespaces
                XCMD(proc_cmd.strip(), wait = WAIT_TILL_DONE) # process
                ### If wanted NON-FT-PROCESSED DATA:
                # XCMD('xtrf n', wait = WAIT_TILL_DONE) # process
                
            # log.info("Finished proc.")
                
def calibrate_1D1H(expno):
    """SREF calibrate 1D1H"""
    XCMD('re '+str(expno), wait = WAIT_TILL_DONE) # read
    XCMD('solvent H2O+D2O_IVTNMR', wait = WAIT_TILL_DONE) # this has a narrower Ref-peak Width range for search
    XCMD('s solvent H2O+D2O_IVTNMR', wait = WAIT_TILL_DONE) # this has a narrower Ref-peak Width range for search
    SREF()
    sr = GETPAR('SR')
    sf = GETPAR('SF')
    log.info('SR / SF after auto-SREF = %s / %s' % (sr, sf))

def get_SR_from_1D1H(expno):
    """Read 1H SR value"""
    XCMD('re '+str(expno), wait = WAIT_TILL_DONE) # read

    # sr = GETPAR('SR') # 20190118. does not work anymore
    sf = GETPAR('SF')
    sr = (float(sf)-float(GETPAR('BF1')))*1e6 # 20190118
    log.info('SR / SF = %s / %s' % (sr, sf))
    return sr,sf

def get_SI_STSI_STSR(expno):
    """Read spectrum size params from specific expno"""
    XCMD('re '+str(expno), wait = WAIT_TILL_DONE) # read
    name, expno, procno, curdir = CURDATA() # get current dataset

    si_stsi_stsr = {
    'SI': (GETPAR("2 SI"), GETPAR("1 SI")),
    'STSI': (GETPAR("2 STSI"), GETPAR("1 STSI")),
    'STSR': (GETPAR("2 STSR"), GETPAR("1 STSR")),
    }
    return si_stsi_stsr

def get_phases_2D(expno):
    """Read phase parameters from specific expno - for both channels"""
    XCMD('re '+str(expno), wait = WAIT_TILL_DONE) # read
    name, expno, procno, curdir = CURDATA() # get current dataset

    # print("2 PHC0 = %.4f" % float(GETPAR("2 PHC0"))
    log.info(expno)
    log.info(GETPAR("2 PHC0"))

    phases = {    
    'PHC0': (float(GETPAR("2 PHC0")), float(GETPAR("1 PHC0"))),
    'PHC1': (float(GETPAR("2 PHC1")), float(GETPAR("1 PHC1"))),
    }
    return phases
    
def get_phase_2D_f2(expno):
    """Only F2 channel!"""
    XCMD('re '+str(expno), wait = WAIT_TILL_DONE) # read
    name, expno, procno, curdir = CURDATA() # get current dataset
    
    phc0 = GETPAR("2 PHC0")
    phc1 = GETPAR("2 PHC1")
    return phc0, phc1    

def get_phase_1D(expno):
    """..."""
    XCMD('re '+str(expno), wait = WAIT_TILL_DONE) # read
    name, expno, procno, curdir = CURDATA() # get current dataset

    phc0 = GETPAR("1 PHC0")
    phc1 = GETPAR("1 PHC1")
    return phc0, phc1

def get_param_1D(expno, param):
    """..."""
    XCMD('re '+str(expno), wait = WAIT_TILL_DONE) # read
    return GETPAR(param)

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
    """Takes system path, returns structure with thousand-based expno numbers"""
    all_dirs = get_immediate_subdirectories( os.path.join(curdir,name) )
    # print all_dirs

    dir_numbers = list( int_filter( all_dirs )) # convert to numbers, dropping non-integers
    # print dir_numbers
    
    expno_sets = []

    for i in range(1,10):
        incr = i*1000
        # print i, incr
        # Selecting not until x999 - cuz there m/b ref expts (e.g. 4998, 4999 for 31P) there
        # Included sorting - otherwise not always can rely on last expno being last.
        expno_sets.append(sorted([e for e in dir_numbers if incr+900 >= e >= incr]))
    # print expno_sets
    return expno_sets
    
def copy_self_to_dataset_dir(curdir,dset_name):
    py_filename = os.path.basename(sys.argv[0])
    py_bak_path = os.path.join(curdir,dset_name,(py_filename+'.copy'))
    log.info('Copying %s to %s ...' % (py_filename,py_bak_path))
    shutil.copy(sys.argv[0], py_bak_path)
    
####################################
## Actual run commands
def main():
    name, expno, procno, curdir = CURDATA() # get current dataset
    # datasetFolder = curdir+"/"+name+"/"    
    
    # ### Replicates for the paper
    # dsets = [
    # '190111_IN70c_R02_co-NUP1_303K_600',
    # # '190110_IN71c_SMN1_co-NUP1_303K_600',
    # # '190111_IN72c_SMN2_co-NUP1_303K_600',
    # ]

    # print curdir
    # print(dsets)
    
    ## Just read current dataset as the one to be processed
    dsets = [name]
    
    # start_log()
    log_file = curdir+'/'+name+'/'+'processing.log'
    print ('======= Log file location =======')
    print log_file
    print ('=================================')
    log.basicConfig(
        filename=log_file,
        # filemode='w', # overwrite file contents instead of appending (default mode is 'a')
        level=log.DEBUG, # minimal level of severity to keep in log: DEBUG INFO WARNING ERROR CRITICAL
        ## if want to debug something specific - use this level + log.warning('...') in the code
        # level=log.WARNING,
        format='%(asctime)s %(levelname)-10s %(message)s',
        #datefmt='%Y%m%d_%H%M%S' # format w/o dashes, but cant easily add millisec
        )

        ## Add log display to the console
    log.getLogger().addHandler(log.StreamHandler())
    
    copy_self_to_dataset_dir(curdir,name) # This contains logging - needs to be run AFTER creation of log object.
    
    log.info('\nThis script allows to process multiple datasets: if "dsets" variable defined in main()')
    log.info('\n\n\n== Starting new logging session')

    n_sets = len(dsets)
    
    for dset in dsets:
        log.info('\n\n== Processing dset %s' % dset)
                
        # Skip datasets which are not sorted or lack 2DHN (4000)
        ### TODO - something weird m/b going on here, cuz it did not work as "not (x or y)" - so had to be implemented as "not x or not y". Later fixed - seems working fine!
        if not (os.path.exists(os.path.join(curdir,dset,'2000')) or os.path.exists(os.path.join(curdir,dset,'2001')) or os.path.exists(os.path.join(curdir,dset,'4000')) or os.path.exists(os.path.join(curdir,dset,'4001'))):
            log.info('Some of expnos used for phase/SR reading (2000,2001,4000,4001) are missing - CHECK IF DSET IS ALREDY SORTED! Skipping this dset for the moment.')
            continue

        read_expt(dset,'1','1',curdir)

        ## Can try to read SR parameter from 1D1H experiment
        if 0: # DO_AUTOMATIC_SREF: not used anymore!
            proc_1D1H_for_SREF(ref_spectrum_for_1D1H_SR)
            calibrate_1D1H(ref_spectrum_for_1D1H_SR)
        
        ## Read the SR and SF parameters
        sr_1H, sf_1H = get_SR_from_1D1H(ref_spectrum_for_1D1H_SR)
                
        log.info('sr_1H=%s sf_1H=%s' % (sr_1H, sf_1H))
        
        ## 20190118 - debugging
        # log.info( type(sr_1H) )
        # EXIT()
        
        if float(sr_1H) == 0.0:
            log.warning('=============================================')
            log.warning('==================== !!! ====================')
            log.warning('Likely need manual SR check in dset %s' % dset)
            log.warning('=============================================')
            log.warning('=============================================')

        expno_sets = get_recorded_expnos(curdir,dset)        
        # print expno_sets

        expnos_2DHN = expno_sets[4-1] # indexing starts with 0
        expnos_31P = expno_sets[5-1] # indexing starts with 0
        expnos_1D1H = expno_sets[2-1] # indexing starts with 0
        expnos_iminos = expno_sets[6-1] # indexing starts with 0        
                                
        # expnos_2DHN = [4000] # TMP - just process 4000 with the phase from 4001
        # print(expno_sets)
        # EXIT()

        ## Processing parameters
        #==============================================
        # Dictionary. Mind that when reading the params F2 will have index "0", F1 - index "1"
        # Each element can contain two entries - can index for multiple dimensions of the spectrum (F1/F2)!!
        # For 1D spectra - only the first entry is used.        
        
        #####################
        ########  31P  ######
        #####################
        if f_process_1D31P and len(expnos_31P) != 0:
            proc_params_31P = {
            'WDW': ('EM',''),
            'LB': (2,''),
            'SI': (65536,''),
            'ABSF1': (14,''),
            'ABSF2': (-26,''),
            'ABSG': (5,''),
            'BC_mod': ('quad',''),
            }
            target_procno='1'
            proc_string_31P = 'efp;apks;absn'
            log.info('\nWill process 31P with following params (expnos, target_procno, SR, proc_params, proc_string):')
            log.info(expnos_31P)
            log.info(target_procno)
            log.info(float(sr_1H)*DSS_SR_1H31P)
            log.info(proc_params_31P)
            log.info(proc_string_31P)
        
            ### Actual processing
            proc_1D(expnos_31P, target_procno, float(sr_1H)*DSS_SR_1H31P, proc_params_31P, proc_string_31P)
        

        ########################
        ########  Iminos  ######
        ########################
        if f_process_iminos and len(expnos_iminos) != 0:
            # get reference phase
            # phc0,phc1 = get_phase_1D(expnos_iminos[0]) # first experiment
            phc0,phc1 = get_phase_1D(expnos_iminos[-1]) # last experiment
            # print(phc0,phc1)
            proc_params_iminos = {
            'PHC0': (phc0, ''),
            'PHC1': (phc1, ''),
            'SI': (32768,''),
            'WDW': ('EM',''),
            'LB': (2,''),
            'BC_mod': ('qpol',''),
            }
            target_procno='1'
            proc_string_iminos = 'fp'
            log.info('\nWill process iminos with following params (expnos, target_procno, SR, proc_params, proc_string):')
            log.info(expnos_iminos)
            log.info(target_procno)
            log.info(float(sr_1H))
            log.info(proc_params_iminos)
            log.info(proc_string_iminos)
        
            ### Actual processing
            proc_1D(expnos_iminos, target_procno, float(sr_1H), proc_params_iminos, proc_string_iminos)

        ########################
        ########  1D1H  ########
        ########################
        if f_process_1D1H and len(expnos_1D1H) != 0:
            # get reference phase
            phc0,phc1 = get_phase_1D(ref_spectrum_for_1D1H_phase)            
            param_bc_mod = get_param_1D(ref_spectrum_for_1D1H_phase,'BC_mod')            
            proc_params_1D1H = {
            'PHC0': (phc0, ''),
            'PHC1': (phc1, ''),
            'SI': (65536,''),
            'WDW': ('EM',''),
            'LB': (2,''),
            'BC_mod': (bc_mod_dict[param_bc_mod],''),
            }
            target_procno='1'
            proc_string_1D1H = 'efp'
            log.info('\nWill process 1D1H with following params (expnos, target_procno, SR, proc_params, proc_string):')
            log.info(expnos_1D1H)
            log.info(target_procno)
            log.info(float(sr_1H))
            log.info(proc_params_1D1H)
            log.info(proc_string_1D1H)
        
            ### Actual processing
            proc_1D(expnos_1D1H, target_procno, float(sr_1H), proc_params_1D1H, proc_string_1D1H)
        
        ########################
        ########  2DHN  ########
        ########################
        if f_process_2DHN and len(expnos_2DHN) != 0:        
            proc_params_2DHN = {
            'WDW': ('QSINE', 'QSINE'),
            'SSB': (2, 2),
            # 'LB': (0, 0),
            # 'FT_mod': ('no', 'no'),
            # 'SI': (4096, 256),
            # 'STSI': (1164, 250),
            # 'STSR': (533, 6),
            # 'PH_mod': ('no', 'no'), 
            'ABSF1': (11, 200),
            'ABSF2': (5.5, 80),
            'ABSG': (5, 5),
            'BC_mod': ('qpol', 'no'),
            }

            phases = get_phases_2D(ref_spectrum_for_2DHN_phase_and_STSI)
            proc_params_2DHN['PHC0'] = phases['PHC0']
            proc_params_2DHN['PHC1'] = phases['PHC1']
            log.info('Phases from ref spectrum %s: phc0(f2/f1)=%s / %s phc1(f2/f1)=%s / %s' % (ref_spectrum_for_2DHN_phase_and_STSI, str(phases['PHC0'][0]), str(phases['PHC0'][1]), str(phases['PHC1'][0]), str(phases['PHC1'][1])))
        
            sss = get_SI_STSI_STSR(ref_spectrum_for_2DHN_phase_and_STSI)
            proc_params_2DHN['SI'] = sss['SI']
            proc_params_2DHN['STSI'] = sss['STSI']
            proc_params_2DHN['STSR'] = sss['STSR']
            log.info('(SI) / (STSI) / (STSR) from ref spectrum: (%s, %s) / (%s, %s) / (%s, %s)' % (str(sss['SI'][0]), str(sss['SI'][1]), str(sss['STSI'][0]), str(sss['STSI'][1]), str(sss['STSR'][0]), str(sss['STSR'][1])))        
        
            target_procno='1'
            proc_string_2DHN = 'xfb n;abs2;abs1'
            log.info('\nWill process 2DHN with following params (expnos, target_procno, SR, proc_params, proc_string):')
            log.info(expnos_2DHN)
            log.info(target_procno)
            log.info(float(sr_1H))
            log.info(proc_params_2DHN)
            log.info(proc_string_2DHN)
        
            ### Actual processing
            proc_2DHN(expnos_2DHN, target_procno, float(sr_1H), proc_params_2DHN, proc_string_2DHN)
        
    log.info('Processing of %i sets completed in %.4f sec.' % (n_sets, time.time()-start_time))
    print ('======= Log file location =======')
    print log_file
    print ('=================================')
    
    # text_editor = '/usr/local/bin/mate'
    # call([text_editor, log_file])    

#####################
# Fore debugging use constructions like:
#print("VARNAME = " + str(VAR))
#print("VARNAME = %.4f" % VAR)
#EXIT() # make a break-point


#####################
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
		main()