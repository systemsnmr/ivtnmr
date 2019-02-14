# Finds and reads expt based on dataset name pattern.
# First argument - name pattern
# Second argument (optional) - target expno

# @Yaroslav Nikolaev. 2017..

import os
import glob

NAME, expno, procno, CURDIR = CURDATA()
#in_expt = 'IN60a'
#expno = 6000
#procno = 1

if len(sys.argv) > 1:
	in_expt = sys.argv[1]

if len(sys.argv) > 2:
	expno = int(sys.argv[2])

if len(sys.argv) > 3:
	procno = int(sys.argv[3])

#print 'in_expt', in_expt

searchpath = os.path.join(CURDIR, ('*%s*' % in_expt) )
#print 'searchpath', searchpath

paths = glob.glob(searchpath)
#print 'paths', paths

#fullpath = [paths[0], str(expno), str(procno), CURDIR]
#print 'fullpath', fullpath

# Cannot use just RE(..) here cuz paths[0] already contains full path!
RE_PATH(os.path.join(paths[0],str(expno),'pdata',str(procno)), show="y") # Syntax/defaults: RE(dataset = None, show = "y")
#'''