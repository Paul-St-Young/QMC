#!/usr/bin/env python
import os
import sys
import subprocess as sp

def backup_target_list(target_dir,name_fmts):
    ''' generate a list of files to backup '''
    targets = []
    for name_fmt in name_fmts:
        proc = sp.Popen('find %s -name "%s" | grep -v %s' % (target_dir,name_fmt,sys.argv[0])
            ,shell=True,stdout=sp.PIPE,stderr=sp.PIPE)
        out,err = proc.communicate()

        targets += out.split('\n')[:-1]
    # end for name_fmt
    return targets
# end def

def backup_in_dir(backup_dir,targets):

    for target in targets:
        subdir = os.path.dirname(target)
        target_dir = backup_dir+'/'+subdir
        if not os.path.exists(target_dir):
            sp.call(['mkdir','-p',target_dir])
        # end if
        failed = sp.call(['cp',target,target_dir])
        if failed:
            raise RunTimeError( "cp failed for %s to %s" % (target,target_dir) )
        # end if
    # end for
    
# end def

if __name__ == '__main__':
    target_dir = '.'
    source_name= sys.argv[1]
    name_fmts  = ['*.py','*.ipnb']

    targets = backup_target_list(target_dir,name_fmts)

    backup_dir = 'scripts_%s' % source_name

    backup_in_dir(backup_dir,targets)
# end __main__
