#!/usr/bin/env python

"""
Miscellaneous system-level utilities
"""
import fnmatch, logging, os, subprocess
from glob import glob
    
logging.basicConfig()

def logged_write(outfilename, mystring, logger = logging, mode = 'w'):
    """
    Writes string to file.

    Includes appropriate logging information.

    @param outfilename Output file name
    @param mystring String to write
    @param logger Logger class
    @param mode Mode to pass to open(). Default: 'w'
    """
    try:
        f = open(outfilename, mode)
        f.write(mystring)
        f.close()
    except IOError, e:
        logger.critical('Could not write to file '+outfilename)
        logger.critical('The data to be written were\n' + mystring)
        raise e
    logger.info('Successfully wrote file '+outfilename)

def expandfilelist(*filelists):
    """
    Returns all files matching a variable number of pathname expressions

    The pathname expression may contain wildcards such as * or ? or others
    parseable by glob.glob(). Returns a unique list of filenames, i.e. files
    matching more than one expression will not be duplicated.

    @param filelists,... Variable number of strings specifying pathname
        expressions

    @returns A sorted list of files that match at least one of the filelists
    """
    #if len(filelists) == 0:
    #   return []
    logger = logging.getLogger('sysutils.expandfilelist')
    logger.debug('Received pathname expressions '+ str(filelists))
    files = set()
    for filelist in filelists:
        try:
            files = files.union(glob(filelist))
        except TypeError:
            logger.warning('Ignoring unparseable glob: '+str(filelist))
    
    files = sorted(files)
    logger.debug('Sorted expanded file list: '+' '.join(files))
    return files

def recursive_glob(pathexpr, root = '.'):
    """
    Recursively globs files matching an expression
    @param pathexpr Path expression to match
    @param root Root to start recursion from. Default: '.'
    @returns list of matching file names
    """
    results = []
    for base, _, files in os.walk(root):
        goodfiles = fnmatch.filter(files, pathexpr)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results

def mostrecent_recursive(pathexpr, root = '.'):
    """
    Recursively finds most recent file matching an expression
    @param pathexpr Path expression to match
    @param root Root to start recursion from. Default: '.'
    @returns name of most recently modified file, else None if no matches
    """
    files = recursive_glob(pathexpr, root)
    if len(files) == 0:
        return
    return max(files, key=os.path.getmtime)

def _exec(command, print_to_screen = False, logfile = None, stdin = None, 
    print_command = True):
    """
    Runs command line using subprocess, optionally returning stdout
    @todo document
    @todo move to sysutils
    """
    print_to_file = (logfile != None)
    if print_to_file:
        f = open(logfile, 'a')
    if print_command:
        print "Executing process: \x1b[1;92m%-50s\x1b[0m Logfile: %s" % (command, logfile)
        if print_to_file:
            print >> f, "Executing process: %s" % command
    if stdin == None:
        p = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        Output, _ = p.communicate()
    else:
        p = subprocess.Popen(command, shell=True, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        Output, _ = p.communicate(stdin)
    if print_to_screen:
        print Output
    if logfile != None:
        f.write(Output)
        f.close()
    return Output

