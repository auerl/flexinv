#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Collection of routines to handle legacy Fortran 77
codes, mostly from Woodhouse and Yu Gu, to setup the
linear system matrices for joint surface and body 
wave tomography at the scale of the entire mantle

:copyright:

Ludwig Auer (ludwig.auer@gmail.com), 2014

:license:

GNU General Public License, Version 3
(http://www.gnu.org/copyleft/gpl.html)


"""


from pylab import *
import ConfigParser
import argparse
import itertools as it
import sys, getopt
import pylab as pl
import numpy as np
import select, re
import time, math
import subprocess
import shutil, os
import glob

# ------------------------------------------------
# subroutine main
# ------------------------------------------------
def main(argv):

    # ------------------------------------------------
    inparam = 'inparam' # Default name of inparam file

    # ------------------------------------------------
    # Get command line options   
    parser = argparse.ArgumentParser()
    parser.add_argument('--inv_dirs',nargs='*',type=str,default=[])
    parser.add_argument('--proj_dirs',nargs='*',type=str,default=[])    
    parser.add_argument('--inv_id', type=str)
    parser.add_argument('--proj_id', type=str)
    parser.add_argument('--rdmp_id', type=str)
    parser.add_argument('--ddmp_id', type=str)
    parser.add_argument('--weight_id', type=str)    
    parser.add_argument('--crust10', action='store_true')
    parser.add_argument('--inversion', action='store_true')
    parser.add_argument('--bowa_matrix', action='store_true')
    parser.add_argument('--suwa_matrix', action='store_true')    
    parser.add_argument('--projection', action='store_true')
    parser.add_argument('--adaptive', action='store_true')
    parser.add_argument('--preani', type=float)
    parser.add_argument('--preani_model', type=str)        
    parser.add_argument('--synth', type=str)    
    parser.add_argument('--xivoigt', action='store_true')
    parser.add_argument('--fmonly', action='store_true')
    parser.add_argument('--pwaves', action='store_true')
    parser.add_argument('--inv_pars', type=str,default=[])
    parser.add_argument('--sw_data_id', type=str)
    parser.add_argument('--bw_data_id', type=str)
    parser.add_argument('--eqincr', type=float)
    parser.add_argument('--layers', type=str)
    parser.add_argument('--qsub', action='store_true')
    parser.add_argument('--bsub', action='store_true')
    parser.add_argument('--ncpu', type=int)
    parser.add_argument('--ophr', type=int)     
    parser.add_argument('--nentries', type=int)
    parser.add_argument('--walltime', type=int) 
    parser.add_argument('--memory', type=int)   
    args = parser.parse_args()
    print ""
    print ""
    print ""
    print "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"    
    print ""
    print " _|_|_|_|  _|        _|_|_|_|  _|      _|  _|_|_|  _|      _|  _|      _|  "
    print " _|        _|        _|          _|  _|      _|    _|_|    _|  _|      _|  "
    print " _|_|_|    _|        _|_|_|        _|        _|    _|  _|  _|  _|      _|  "
    print " _|        _|        _|          _|  _|      _|    _|    _|_|    _|  _|    "
    print " _|        _|_|_|_|  _|_|_|_|  _|      _|  _|_|_|  _|      _|      _|      "
    print ""
    print ""
    print "- - - - - Welcome to setup.py: - - - - - - "
    print ""
    print ""
    if args.suwa_matrix:
        print "----------------------------------"        
        print "setup.py: surface wave matrix mode"
        print "----------------------------------"
    if args.bowa_matrix:
        print "----------------------------------"                
        print "setup.py: body wave matrix mode"
        print "----------------------------------"                
    if args.projection:
        print "----------------------------------"                
        print "setup.py: projection mode"
        print "----------------------------------"                
    if args.inversion:        
        print "----------------------------------"                
        print "setup.py: inversion mode"
        print "----------------------------------"                


    # ------------------------------------------------
    # Computational parameters related to the executing machine
    machine={}    
    machine["sche"] = readip(inparam,'machine','sche','list')[0] # scheduler
    machine["ophr"] = readip(inparam,'machine','ophr','int') # data per hr (estimate)
    machine["ncpu"] = readip(inparam,'machine','ncpu','int') # nrs of cpus    
    machine["nent"] = readip(inparam,'machine','nent','int')
    machine["wall"] = readip(inparam,'machine','wall','int')
    machine["memo"] = readip(inparam,'machine','memo','int')    
    
    if args.ophr: machine["ophr"]=args.ophr
    if args.nentries: machine["nent"]=args.nentries
    if args.ncpu: machine["ncpu"]=args.ncpu
    if args.memory: machine["memo"]=args.memory
    if args.walltime: machine["wall"]=args.walltime
    if args.bsub: machine["sche"]='bsub'        
    if args.qsub: machine["sche"]='qsub'        
    if args.bsub and args.qsub:
        print "Both lsf and torque specified, what do you want?"
        print "Preparing submit scripts for a workstation!"
        machine["sche"] ='pc'


    # ------------------------------------------------
    # Read basic matrix parameters
    eqincr  = readip(inparam,'param','eqincr','float')
    xivoigt = readip(inparam,'param','xivoigt','bool')    
    adaptive= readip(inparam,'param','adaptive','bool')

    pwaves=False
    if args.pwaves: pwaves=True
    if args.eqincr: eqincr=args.eqincr
    if args.adaptive: adaptive=True
    if args.xivoigt: xivoigt=True
    if adaptive:
        if eqincr!=1.25:
            print 'Inconsistency: Adaptive, but eqincr'
            print 'is',eqincr,'-> setting it to 1.25'
            eqincr=1.25


    # ------------------------------------------------
    # Figure out which scheduler is used

    # ------------------------------------------------
    # Read basic folder parameters
    rootdir = readip(inparam,'folders','root','list')[0]
    kerdir  = readip(inparam,'folders','kernel','list')[0]
    matdir  = readip(inparam,'folders','matrix','list')[0]
    outdir  = readip(inparam,'folders','output','list')[0]        

    # ------------------------------------------------
    # Read which data to process
    swdatset = readip(inparam,'swdata','set','list')[0]
    bwdatset = readip(inparam,'bwdata','set','list')[0]
    if args.sw_data_id: swdatset=args.sw_data_id
    if args.bw_data_id: bwdatset=args.bw_data_id
    swdatdir = readip(inparam,'swdata',swdatset,'list')[0]
    bwdatdir = readip(inparam,'bwdata',bwdatset,'list')[0]   

    # ------------------------------------------------
    # Read radial parameterization
    layset = readip(inparam,'layers','set','list')[0]
    if args.layers: layset=args.layers
    layers = readip(inparam,'layers',layset,'list')
    nlay   = len(layers)
          
    # ------------------------------------------------
    # Read weighting and damping schedules
    pharr,wharr = read_weighting_scheme(inparam,args.weight_id)
    rdamp = readip(inparam,'damping','rdamp','list')
    ndamp = readip(inparam,'damping','ndamp','list')
    ddamp = readip(inparam,'damping','ddamp','list')
    
    # ------------------------------------------------
    # Read basic inversion parameters
    synth = readip(inparam,'invparam','synth','bool')
    model = readip(inparam,'invparam','model','list')[0]
    if args.synth:
        synth = True
        model = args.synth
    cfacsw= readip(inparam,'invparam','cutoffsw','float')
    cfacbw= readip(inparam,'invparam','cutoffbw','float')    

    # ------------------------------------------------
    # Prescribed anisotrop
    preani=0.
    preani_model=" "
    if args.preani:
        preani=args.preani
        preani_model=args.preani_model

    # ------------------------------------------------
    # Read which data to combine in inversion mode

    if args.inv_dirs:
        invdirs = args.inv_dirs
    else:
        try:
            invdirset = readip(inparam,'invparam','set','list')[0]
            invdirs = readip(inparam,'invparam',invdirset,'list')
        except:
            print "ERROR:"
            print "You have not provided matrix folders via the command line"
            print "and I could not read them from the inparam file... "
            sys.exit()
        
    invdirs = np.array(invdirs).reshape((len(invdirs)/3,3))    
    bwupscale=[]; bwfolders=[]
    swupscale=[]; swfolders=[]
    adpupscale=[]; adpfolders=[]            
    for i in range(0,len(invdirs)):
        if invdirs[i,0]=='B':
            bwfolders.append(invdirs[i,1])
            bwupscale.append(invdirs[i,2])
        elif invdirs[i,0]=='S':
            swfolders.append(invdirs[i,1])
            swupscale.append(invdirs[i,2])
        elif invdirs[i,0]=='A':
            adpfolders.append(invdirs[i,1])
            adpupscale.append(invdirs[i,2])            
        else:
            print 'some error'
            sys.exit()

    # -----------------------------------------------
    # Read physical inversion parameters
    invpars = []
    invptmp = readip(inparam,'param','invpars','list')
    invptmp = [(int(el.strip('[').strip(',').strip(']')))
                for el in invptmp]    
    if args.inv_pars:
        invptmp = np.array(args.inv_pars.split())
    for i in invptmp:
        invpars.append(int(i))

    
    # ------------------------------------------------
    # Read projection parameters
    olddws  = readip(inparam,'projparam','olddws','bool')
    dwsfile = readip(inparam,'projparam','dwsfile','list')[0]
    thresh  = readip(inparam,'projparam','thresh','list')    

    # ------------------------------------------------
    # Read which data to combine in projection mode
    projdirset = readip(inparam,'projparam','set','list')[0]    
    projdirs = readip(inparam,'projparam',projdirset,'list')
    if args.proj_id and not args.proj_dirs:
        projdirset=args.proj_id
        projdirs = readip(inparam,'projparam',projdirset,'list')        
    elif args.proj_id and args.proj_dirs:
        projdirs = args.proj_dirs
        projdirset=args.proj_id
    elif not args.proj_id and args.proj_dirs:              
        projdirs=args.proj_dirs
        projdirset=raw_input('Enter project id for \
                              this projection:')
    projdirs = np.array(projdirs).reshape((len(projdirs)/3,3))
    bwprojupscale=[]; bwprojfolders=[]
    swprojupscale=[]; swprojfolders=[]        
    for i in range(0,len(projdirs)):
        if projdirs[i,0]=='B':
            bwprojfolders.append(projdirs[i,1])
            bwprojupscale.append(projdirs[i,2])
        elif projdirs[i,0]=='S':
            swprojfolders.append(projdirs[i,1])
            swprojupscale.append(projdirs[i,2])                
        else:
            print 'some error'
            sys.exit()
            
    # ------------------------------------------------
    # Find names of kernel/matrix output folders
    ipar='-p'    
    for i in invpars:
        ipar=ipar+str(i)
    if xivoigt:  # This part is still a separate fork
        ipar='234'
        av='-av'
    else:
        av=''
        
    # ------------------------------------------------
    # Set up kernel and matrix directory paths
    if args.crust10:        
        keroutdir = kerdir+'/c1'+av+'-'+str(nlay)+'/'
    else: # use crust 2.0
        keroutdir = kerdir+'/c2'+av+'-'+str(nlay)+'/'
        
    bwoutdir  = matdir+'/'+bwdatset+av+ipar+'-'+ \
        str(nlay)+'-'+str(eqincr)+'/'
    swoutdir  = matdir+'/'+swdatset+av+ipar+'-'+ \
        str(nlay)+'-'+str(eqincr)+'/'
    
    keroutdir = keroutdir.replace('//','/')
    swoutdir = swoutdir.replace('//','/')
    bwoutdir = bwoutdir.replace('//','/')


    # ------------------------------------------------
    # Check if ker/mat/inv already exist and run the subroutines
    if args.suwa_matrix:
        
        # Does surface wave kernel folder exist?
        iswkers=check_dir(keroutdir) # -> True if it iexists
        if not iswkers: # Otherwise (iskers=false), delete?
            yesno=raw_input('... delete? (yes/no)')
            if yesno=='yes' or yesno=='YES': # if 'n' iswkers=False
                os.system('rm -rf '+keroutdir)
                make_dir(keroutdir)
                iswkers=True
        if iswkers: # Compute them
            print 'Preparing surface wave kernel local profiles'
            print 'This may take a couple of minutes!'            
            if args.crust10:
                prepare_suwa_kernels_crust10(keroutdir,rootdir,
                                         xivoigt,layers,machine)
            else:
                prepare_suwa_kernels_crust20(keroutdir,rootdir,
                                         xivoigt,layers,machine)
            

        # Does surface matrix folder exist?
        iswmatr=check_dir(swoutdir) # -> gets True if it doesn't exist and gets created
        if not iswmatr: # Otherwise (iswmatr=False) it exists, so ask whether delete...
            yesno=raw_input('... delete? (y/n)')
            if yesno=='y' or yesno=='Y': # if 'n' iswmatr will stay False
                os.system('rm -rf '+swoutdir)
                make_dir(swoutdir)
                iswmatr=True                
        if iswmatr: # Compute them!
            print 'Preparing submit script for sw matrices!'
            make_dir(swoutdir+'dws/')
            prepare_suwa_matrix(swdatdir,swoutdir,keroutdir,
                                rootdir,eqincr,nlay,adaptive,
                                swdatset,args.fmonly,invpars,
                                xivoigt,machine)
            
            # Save important parameters
            params = ConfigParser.RawConfigParser()
            params.add_section('Parameters')
            params.set('Parameters','dataid',swdatset)
            params.set('Parameters','eqincr',eqincr)
            params.set('Parameters','adaptive',adaptive)
            params.set('Parameters','xivoigt',xivoigt)
            params.set('Parameters','layers',layers)
            params.set('Parameters','invpars',invpars)            
            paramfile=open(swoutdir+'/param', 'wb')
            params.write(paramfile)
            suwa = open(matdir+'/suwa','a+')
            suwa.write(swoutdir+'\n')
            suwa.close()
                            
    # ------------------------------------------------
    if args.bowa_matrix:
        ibwmatr = check_dir(bwoutdir)
        if not check_dir(bwoutdir):
            yesno=raw_input('Overwrite? (y/n)')
            if yesno=='y' or yesno=='Y':
                os.system('rm -rf '+bwoutdir)
                make_dir(bwoutdir)
                ibwmatr = True
        if ibwmatr:
            print 'Computing body wave matrices!'
            make_dir(bwoutdir+'dws/')
            make_dir(bwoutdir+'tmp/')            
            # Save important parameters
            params = ConfigParser.RawConfigParser()
            params.add_section('Parameters')
            params.set('Parameters','dataid',bwdatset)
            params.set('Parameters','eqincr',eqincr)
            params.set('Parameters','adaptive',adaptive)
            params.set('Parameters','xivoigt',xivoigt)
            params.set('Parameters','layers',layers)
            params.set('Parameters','invpars',invpars)            
            paramfile=open(bwoutdir+'/param', 'wb')
            params.write(paramfile)
            prepare_bowa_matrix(bwdatdir,bwoutdir,rootdir,eqincr,
                                bwdatset,layers,adaptive,xivoigt,
                                invpars)                    
            bowa = open(matdir+'/bowa','a+')
            bowa.write(bwoutdir+'\n')
            bowa.close()
                    
    # ------------------------------------------------            
    if args.projection:

        # Read parameters from inputfolders
        aav=''
        allfolders=[]
        allfolders.append(bwprojfolders)
        allfolders.append(swprojfolders)
        allfolders=[item for sublist in allfolders for item in sublist]
        print "Processing folders:", allfolders        
        projnlay = len(readip(allfolders[0]+'/param',
                     'Parameters','layers','list'))
        projav = readip(allfolders[0]+'/param',
                     'Parameters','xivoigt','bool')
        projip = readip(allfolders[0]+'/param',
                     'Parameters','invpars','list')
        projip=[(int(el.strip('[').strip(',').strip(']')))
               for el in projip]
        projad = readip(allfolders[0]+'/param',
                     'Parameters','adaptive','bool')

        # Perform a number of sanity checks
        params = ConfigParser.ConfigParser(); 
        for para in ['layers','eqincr','xivoigt','invpars','adaptive']:
            check=[]            
            for item in allfolders:
                check.append(readip(item+'/param','Parameters',para,'list'))
            if not all_same(check):
                print 'Parameter inconsistency at parameter', para
                sys.exit()

        # Some more sanity checks
        if projnlay != nlay:
            print 'Parameter inconsistency at parameter layers', layers
            sys.exit()
        if projav != xivoigt:
            print 'Parameter inconsistency at parameter xivoigt', xivoigt
            sys.exit()
        if projad != adaptive:
            print 'Parameter inconsistency at parameter adaptive', adaptive
            sys.exit()
        if projip != invpars:
            print 'Parameter inconsistency at parameter invpars', invpars
            sys.exit()            
        if not projad:
            print 'Parameter inconsistency at parameter projad', projad
            sys.exit()
        if projav: aav='-av'

        # Name project and store parameter files
        projident = raw_input('Enter an ID to dub the projection (enter for default): ')        
        if projident=='':
            projoutdir  = matdir+'/'+projdirset+'-adp'+aav+'-'+ \
                str(nlay)+'-'+str(eqincr)+'/'         
        else:
            projoutdir  = matdir+'/'+projident+'-adp'+aav+'-'+ \
                str(nlay)+'-'+str(eqincr)+'/'
        projoutdir = projoutdir.replace('//','/')                            
        doproj=check_dir(projoutdir)        
        if not doproj:
            yesno=raw_input('Overwrite? (y/n)')
            if yesno=='y' or yesno=='Y':
                os.system('rm -rf '+projoutdir)
                make_dir(projoutdir)

        # Finally, create projection submit script
        prepare_projection(swprojfolders,bwprojfolders,
                           projoutdir,pharr,wharr,
                           swprojupscale,bwprojupscale,
                           synth,model,cfacsw,
                           cfacbw,rdamp,ddamp,
                           ndamp,nlay,eqincr,
                           rootdir,thresh,olddws,
                           dwsfile,nentries,
                           invpars,ipar,machine)
        

    # ------------------------------------------------
    if args.inversion:
        # Compute synthetics if requested


        # Read parameters from inputfolders        
        iad=''; iav=''
        allfolders=[]
        allfolders.append(bwfolders)
        allfolders.append(swfolders)
        allfolders.append(adpfolders)                
        allfolders=[item for sublist in allfolders for item in sublist]
        print "Processing folders:", allfolders        
        invnlay = len(readip(allfolders[0]+'/param',
                     'Parameters','layers','list'))
        invav = readip(allfolders[0]+'/param',
                     'Parameters','xivoigt','bool')
        invip = readip(allfolders[0]+'/param',
                     'Parameters','invpars','list')        
        invip=[(int(el.strip('[').strip(',').strip(']')))
               for el in invip]
        invad = readip(allfolders[0]+'/param',
                     'Parameters','adaptive','bool')


        # Perform a number of consistency checks
        params = ConfigParser.ConfigParser(); 
        for para in ['layers','eqincr','xivoigt','invpars','adaptive']:
            check=[]            
            for item in allfolders:
                check.append(readip(item+'/param','Parameters',para,'list'))
            if not all_same(check):
                print 'Parameter inconsistency at parameter', para
                sys.exit()

        # Some more consistency checks
        if invnlay != nlay:
            print 'Parameter inconsistency at parameter layers', layers
            sys.exit()
        if invav != xivoigt:
            print 'Parameter inconsistency at parameter xivoigt', xivoigt
            sys.exit()
        if invip != invpars:
            print 'Parameter inconsistency at parameter invpars', invpars
            sys.exit()            
        if invad != adaptive:
            print 'Parameter inconsistency at parameter adaptive', adaptive
            sys.exit()                                   
        if invad: iad='adp.'
        if invav: iav='av.'

        # Read damping schedule
        dmpscheme_r = read_damping_scheme(inparam,invnlay,args.rdmp_id,5)
        dmpscheme_d = read_damping_scheme(inparam,invnlay,args.ddmp_id,2)
        
        # Name project and store parameter files
        invident=''
        if not args.inv_id:
            invident = raw_input('Enter a project ID for this inversion (enter for default):')
        if args.inv_id:
            invident = args.inv_id


        rdid=readip(inparam,'damping','rdamp_scheme','list')[0]
        ddid=readip(inparam,'damping','ddamp_scheme','list')[0]        
        wid=readip(inparam,'weighting','set','list')[0]

        if args.rdmp_id: rdid=args.rdmp_id
        if args.ddmp_id: ddid=args.ddmp_id
        if args.weight_id: wid=args.weight_id

        if invident=='':
            invoutdir = outdir+'/inv.'+rdid+'.'+ddid+'.'\
                +wid+'.'+iav+iad+ipar+'l'+str(nlay)+'.e'\
                +str(eqincr).replace('.0','')+'.r'\
                +str(rdamp[0])+'-'+str(rdamp[len(rdamp)-1])+'/'
        else:
            invoutdir = outdir+'/'+invident+'/'

        invoutdir = invoutdir.replace('//','/')                            
        doinv = check_dir(invoutdir)
        if not doinv:
            yesno=raw_input('Overwrite? (y/n)')
            if yesno=='y' or yesno=='Y':
                os.system('rm -rf '+invoutdir)
                make_dir(invoutdir)
                doinv = True

        # Dump definitions associated with damping
        dmp_file_r = invoutdir+'/rdamping.in' # ASCII to be read by flexinv90
        dmp_file_d = invoutdir+'/ddamping.in' # ASCII to be read by flexinv90
        dump_damping(dmp_file_r,rdamp,ddamp,ndamp,dmpscheme_r,5) # has 5 columns
        dump_damping(dmp_file_d,rdamp,ddamp,ndamp,dmpscheme_d,2) # has 2 columns
        layfil = invoutdir+'layers.in'
        layfil = layfil.replace('//','/')
        dump_layers(layers,layfil)

        # Dump inversion parameters to config file
        params = ConfigParser.RawConfigParser()
        params.add_section('Parameters')
        params.set('Parameters','eqincr',eqincr)
        params.set('Parameters','adaptive',adaptive)
        params.set('Parameters','xivoigt',xivoigt)
        params.set('Parameters','layers',' '.join(str(x) for x in layers))
        params.set('Parameters','rdamp',' '.join(str(x) for x in rdamp))
        params.set('Parameters','ndamp',' '.join(str(x) for x in ndamp))
        params.set('Parameters','ddamp',' '.join(str(x) for x in ddamp))  
        paramfile=open(invoutdir+'/param', 'wb')
        params.write(paramfile)

        #  calls subroutines that prepare the inversions
        if doinv and not adaptive:
            os.system('cp '+inparam+' '+invoutdir)
            prepare_inversion_reg(swfolders,bwfolders,
                                  invoutdir,pharr,wharr,
                                  swupscale,bwupscale,
                                  synth,model,cfacsw,
                                  cfacbw,rdamp,ddamp,
                                  ndamp,nlay,eqincr,
                                  dmp_file_r,dmp_file_d,
                                  rootdir,invpars,preani,
                                  preani_model,machine)
            
        # Inversion of variable grid matrices            
        elif doinv and adaptive:
            os.system('cp '+inparam+' '+invoutdir)
            prepare_inversion_adp(adpfolders,invoutdir,
                                  synth,model,cfacsw,
                                  cfacbw,rdamp,ddamp,
                                  ndamp,nlay,eqincr,
                                  dmp_file_r,dmp_file_d,
                                  rootdir,invpars,machine)
                         


#
#  ---> Below you find all subroutines
#            

        
# ================================================================
#  Subroutines that read/write damping schemes

# Read difference/roughness damping scheme
def read_damping_scheme(fname,nlay,dmpid,ncol):
    dmpset = readip(fname,'damping','rdamp_scheme','list')[0]
    if dmpid: dmpset=dmpid
    dmpsch = readip(fname,'damping',dmpset,'list')
    if (len(dmpsch)/ncol)-1 != nlay:
        print 'Error: nlay and roughness damping scheme inconsistent'
        print len(dmpsch)/ncol-1, nlay
        sys.exit()        
    dmpsch = array(dmpsch).reshape(((nlay+1),ncol))
    return dmpsch

# Dump difference/roughness damping scheme
def dump_damping(fname,rdamp,ddamp,ndamp,dmpsch,ncol):
    dfil=open(fname,'w+')
    [dfil.write(str(el)+' ') for el in rdamp]; dfil.write('\n')
    [dfil.write(str(el)+' ') for el in ndamp]; dfil.write('\n')
    [dfil.write(str(el)+' ') for el in ddamp]; dfil.write('\n')
    dfil.write('#######\n')
    for i in range(0,np.shape(dmpsch)[0]):
        if i==1: dfil.write('#######\n')
        for j in range(0,ncol):
            dfil.write(str(dmpsch[i,j])+' ')
        dfil.write('\n')
    dfil.close()
    

# ================================================================
#  Subroutines that reads weighting scheme
def read_weighting_scheme(fname,weightid):
    pharr = []; wharr = []
    setss = readip(fname,'weighting','set','list')[0]
    if weightid: setss=weightid
    llist = readip(fname,'weighting',setss,'list')            
    for item in llist:
        a = readip(fname,'weighting',item,'list')
        b = len(a)
        if (b % 2) != 0:
            print 'There seems to be an error \
            with your weighting schedule!'
            sys.exit()
        for c in range(0,(b/2)):
            if not 'B' in item:
                phaid = '.'+item+'.'+a[c]+'.'
            else:
                phaid = '.'+a[c]+'.'                    
            pharr.append(phaid)                
            wharr.append(a[c+(b/2)])
    # remove 0 weighted entries
    for item in reversed(wharr):
        if item[-5:] == '0':
            pharr.pop(wharr.index(item))
            wharr.remove(item)
    return pharr,wharr


# ================================================================
#  Subroutines to dump files with layers
def dump_layers(layers,fname):    
    lays=open(fname,'w+')
    lays.write(str(len(layers))+'\n')
    for layer in layers:
        lays.write(str(float(layer))+'\n')
    lays.close()



# ================================================================
#  Subroutines that reads/writes different datatypes
#  from inparam files
def readip(fname,grp,sub,typ):
    config = ConfigParser.ConfigParser()
    config.read(fname)
    if typ == 'list':
        zurueck = config.get(grp,sub).split()
    elif typ == 'int':
        zurueck = config.getint(grp,sub)
    elif typ == 'bool':
        zurueck = config.getboolean(grp,sub)
    elif typ == 'float':
        zurueck = config.getfloat(grp,sub)
    else:
        print 'Strange type, read_inparam'
        sys.exit()        
    return zurueck


# ================================================================
#  Subroutines related to list and string handling
#  Returns true if all items of list are the same
def all_same(items):
    it = iter(items)
    for first in it:
        break
    else:
        return True  # empty case, note all([]) == True
    return all(x == first for x in it)

# pattern search is to put paths in their correct order
def first_substr(thelist, substring):
    for i, s in enumerate(thelist):
        if substring in s:
            return i
    return -1

# divides lists in chungs of n items
def sep_chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]




# ================================================================
#  Subroutines performing basic system calls

# Count lines    
def line_count(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# Create directories
def make_dir(newdir):
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            make_dir(head)
        if tail:
            os.mkdir(newdir)

# check if kernel folder exists, if not create one
def check_dir(newdir):
    if os.path.isdir(newdir):
        print newdir, 'EXISTS already!'
        runcomp = False
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        print newdir, 'Does NOT exist, working ...'
        runcomp = True
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            make_dir(head)
        if tail:
            os.mkdir(newdir)
    return runcomp
            


# ================================================================
#  Subroutine returning nominal nr of inversion parameters
#  mainly used for subsequent sanity checks
        
def get_nombloc(eqincr):
    if eqincr==1.25:
        nombloc=31808        
    elif eqincr==2:
        nombloc=12436
    elif eqincr==3:
        nombloc=5536
    elif eqincr==5:
        nombloc=1988
    elif eqincr==6:
        nombloc=1388
    else:
        print 'Invalid equincr ', eqincr
        exit()
    return nombloc



# ================================================================
#  Major subroutine that prepares submit scripts for surface
#  wave kernel computations. The resulting scripts are meant
#  to be run in parallel on a multi-core machine
#  
def prepare_suwa_kernels_crust10(keroutdir,rootdir,xivoigt,
                                 layers,machine):

    # iterative mode not really used
    itera=0; itera3="%03d" % itera; ixv='n'
    if xivoigt: ixv='y'       
    ncpu=machine["ncpu"]
    
    # set up all 64800 crustal locations
    locations=[]
    for lat in xrange(1,181):
        for lon in xrange(1,361):
            pn1="%03d" % lat; pn2="%03d" % lon
            pixnam=pn1+pn2; locations.append(pixnam)

    # crustal and reference model paths
    kerbindir= rootdir+'/mat/suwa/kernels/'
    refmfile = kerbindir+"/mod/anipremsw"
    cdirname = kerbindir+"/mod/crust10/"
    refmfile = refmfile.replace('//','/')
    cdirname = cdirname.replace('//','/')

    # radial parameterization
    layfil = keroutdir+'layers.in'
    layfil = layfil.replace('//','/')
    dump_layers(layers,layfil)

    # start the crlymod part
    logfil=keroutdir+'crlyrmod.log'; model=[]    
    with open(logfil,"a") as out:
        
        cmdstr=[kerbindir+"bin/crlyrmod_10",
                "-R",refmfile,
                "-C",cdirname,
                "-O",keroutdir]

        proc=subprocess.Popen(cmdstr,stdout=out)            
        proc.wait()

    # start vox2pert to get from .pcn to .pert
    cpu=0; cpu3="%03d" % cpu; count=0
    punchcardname=keroutdir+'ker.'+str(cpu3)+'.sh'   
    punchcard=open(punchcardname,'w+')        
    submitname=keroutdir+'submit.sh'   
    submit=open(submitname,'a')

    # loop over all locations
    for i in range(0,len(locations)):
        if i%648==0:
            print 'Setting up bash scripts: ',(i/648)+1,' %'

        inpuname=keroutdir+'L_'+locations[i]+'_'
        outpname=inpuname+itera3
        if itera==0:
            pertname=cdirname+"L_000000_000.pert"
        else:
            pertname=inpuname+itera3+'.pert'

        if len(inpuname)>80 or len(pertname)>80 \
            or len(outpname)>80 or len(layfil)>80:
            print 'Paths consist of more than 80 chars,\
                   please shorten them. Sorry, this is \
                   a very old software.'
            sys.exit()
        
        cmdstr1=kerbindir+"bin/ptlyrmod -I "+inpuname+" -P "+\
                pertname+" -O "+outpname+" -L "+layfil+"\n"             
        cmdstr2=kerbindir+"bin/min_kernel "+outpname+"\n"
        cmdstr3=kerbindir+"bin/integrate_modes -I "+outpname+\
                " -L "+layfil+" -P "+ixv+"\n"
                
        punchcard.write(cmdstr1)
        punchcard.write(cmdstr2)
        punchcard.write(cmdstr3)
        
        count=count+1
        if count>=(int(64800/ncpu)): #64800 is number of crust 1.0 profiles
            submit.write(punchcardname+' &\n')
            cpu=cpu+1; cpu3="%03d" % cpu; count=0
            punchcard.close() # open punchcard for next cpu     
            punchcardname=keroutdir+'ker.'+str(cpu3)+'.sh'
            punchcard=open(punchcardname,'w+')
            
    punchcard.close()    
    clean=open(keroutdir+'clean.sh','w+')
    cmdstr='rm -rf '+keroutdir+'*.bin *.out ker.*.sh L_*_ L_*_000'
    clean.write(cmdstr)        
    clean.close(); submit.close()
    os.system('chmod 755 '+keroutdir+'/clean.sh '
              +keroutdir+'submit.sh')




# ================================================================
#  Major subroutine that prepares submit scripts for surface
#  wave kernel computations. The resulting scripts are meant
#  to be run in parallel on a multi-core machine
#  

def prepare_suwa_kernels_crust20(keroutdir,rootdir,xivoigt,
                                 layers,machine):    

    # Iterative mode not used, so far
    itera=0; itera3="%03d" % itera; ixv='n'
    if xivoigt: ixv='y'
    ncpu=machine["ncpu"]
    time=machine["wall"]

    # Paths to crustal and reference models
    kerbindir= rootdir+'/mat/suwa/kernels/'        
    swrefmod = kerbindir+'/mod/anipremsw'  # ref mod in a suitable format
    crustdir = kerbindir+'/mod/crust20/'   # folder containing crust 2.0
    regiofil = crustdir+'/Regions_2.0'     # crust 2.0 regions definition

    regiofil = regiofil.replace('//','/')
    crustdir = crustdir.replace('//','/')    
    swrefmod = swrefmod.replace('//','/')

    # Radial parameterization    
    layfil = keroutdir+'layers.in'
    layfil = layfil.replace('//','/')
    dump_layers(layers,layfil)

    # create model from CRUST2.0 and ANIPREM for each loc
    regions=np.genfromtxt(regiofil,dtype=None)

    # ptlyrmod using a reference of output of vox2pert
    reg=[]; loc=[]

    # set up list of 16200 locs and ~360 regs
    for liste in regions:         
        reg.append(liste[4])
        if liste[1]<0:
            lon=(360+int(liste[1]))/int(liste[3])
        else:
            lon=int(liste[1])/int(liste[3])
        pn1num=((90-int(liste[0]))/int(liste[2]))+1
        pn2num=(lon+1)
        pn1="%03d" % pn1num; pn2="%03d" % pn2num
        pixnam=pn1+pn2; loc.append(pixnam)

    # start the crlymod part
    logfil=keroutdir+'crlyrmod.log'; model=[]    
    with open(logfil,"a") as out:
        # loop over regions
        for liste in regions:
            model.append(liste[4])
            model=list(set(model)) # remove dublicates
        # submit many jobs    
        for regio in model:

            if len(swrefmod)>80 or len(crustdir)>80:
                print 'Paths consist of more than 80 chars,\
                       please shorten them. Sorry, this is \
                       a very old software.'
                sys.exit()

            cmdstr=[kerbindir+"bin/crlyrmod_20",
                    "-a","-R",swrefmod,"-C",regio,
                    "-M",crustdir]
            proc=subprocess.Popen(cmdstr,stdout=out)            
            proc.wait()
            if os.path.exists(keroutdir+regio):
                os.remove(keroutdir+regio)
            shutil.move(regio,keroutdir)
            
    # start vox2pert to get from .pcn to .pert
    cpu=0; cpu3="%03d" % cpu; count=0
    punchcardname=keroutdir+'ker.'+str(cpu3)+'.sh'   
    punchcard=open(punchcardname,'w+')        
    submitname=keroutdir+'submit.sh'   
    submit=open(submitname,'a')

    # loop over all locations
    for i in range(0,len(loc)):
        if i%162==0:
            print 'Setting up bash scripts: ',(i/162)+1,' %'
        inpuname=keroutdir+'L_'+loc[i]+'_'
        outpname=inpuname+itera3
        if itera==0:
            pertname=crustdir+"L_000000_000.pert"
        else:
            pertname=inpuname+itera3+'.pert'
        # copy regionfiles to their loc fil
        shutil.copy(keroutdir+reg[i],inpuname)

        if len(inpuname)>80 or len(pertname)>80 \
            or len(outpname)>80 or len(layfil)>80:
            print 'Paths consist of more than 80 chars,\
                   please shorten them. Sorry, this is \
                   a very old software.'
            sys.exit()
        
        cmdstr1=kerbindir+"bin/ptlyrmod -I "+inpuname+" -P "+\
                pertname+" -O "+outpname+" -L "+layfil+"\n"             
        cmdstr2=kerbindir+"bin/min_kernel "+outpname+"\n"
        cmdstr3=kerbindir+"bin/integrate_modes -I "+outpname+\
                " -L "+layfil+" -P "+ixv+"\n"
                
        punchcard.write(cmdstr1)
        punchcard.write(cmdstr2)
        punchcard.write(cmdstr3)
        
        count=count+1
        if count>=(int(16200/ncpu)): #16200 is number of crust 2.0 profiles
            if machine["sche"]=='bsub':
                submit.write('bsub -R lustre -W'+str(time)+':59 < '
                             +punchcardname+'\n')
            elif machine["sche"]=='qsub':
                submit.write('qsub -l nodes=1:ppn=1,walltime='+str(time)
                             +':00:00 '+punchcardname+'\n')
            elif machine["sche"]=='pc':
                submit.write(punchcardname+' &\n')
                            
            cpu=cpu+1; cpu3="%03d" % cpu; count=0            
            punchcard.close() # open punchcard for next cpu     
            punchcardname=keroutdir+'ker.'+str(cpu3)+'.sh'
            punchcard=open(punchcardname,'w+')
            
    punchcard.close()    
    clean=open(keroutdir+'clean.sh','w+')
    cmdstr='rm -rf '+keroutdir+'*.bin *.out ker.*.sh L_*_ L_*_000'
    clean.write(cmdstr)        
    clean.close(); submit.close()
    os.system('chmod 755 '+keroutdir+'/clean.sh '
              +keroutdir+'submit.sh')



# ================================================================
#  Major subroutine that prepares submit scripts for surface
#  wave matrix computations. The resulting scripts are meant
#  to be run in parallel on a multi-core machine
#  

def prepare_suwa_matrix(swdatdir,swoutdir,keroutdir,
                        rootdir,eqincr,nlay,adaptive,
                        swdatset,fmonly,invpars,
                        xivoigt,machine):

    # Some basic parameters 
    swdatun = 1 # Unit of data, 1=phase delay [s], 2=phase velocity
    itera   = 0 # Iteration number, always zereo for the moment
    irorrl  = 1 # rayleigh and love, or only one of them
    ophr    = machine["ophr"]


    cpu=0 # Index for the nr of cpus 
    cpu3="%03d" % cpu

    # vph vpv vsh vsv part is kept general now
    matmode=str(len(invpars))+'\n'
    for i in invpars:
        matmode=matmode+str(i)+'\n'       
    if xivoigt:  # This part is still a separate fork
        matmode='\n2\n3\n4\n'

    # Not very elegant
    if adaptive:
        iadapt=1
    else:
        iadapt=0        
    if fmonly:
        add='_fm'
    else:
        add=''
        
    obsfilsw = swdatdir+'/observed_modes_'
    swbindir = rootdir+'/mat/suwa/'
    obsfilsw = obsfilsw.replace('//','/')
    swdatdir = swdatdir.replace('//','/')
    swbindir = swbindir.replace('//','/')    
    swoutdir = swoutdir.replace('//','/')    
    nombloc = get_nombloc(eqincr)
    
    for typ in ['L','R']:       
        irorrl=irorrl+1 # irorrl: read RL (1), L(2), R(3)
        obslines=open(obsfilsw+typ+add+'.txt')
        obslines=obslines.readlines()
        submitname=swoutdir+'submit.sh'   
        submit=open(submitname,'a')
        for obs in obslines:
            tmpnam=swoutdir+'obs.'+str(cpu3)+'.'
            cputmp=open(tmpnam+typ+'.txt','w+')
            cputmp.write(obs)
            cputmp.close()
            logfil=swoutdir+'sw.'+str(cpu3)+'.log'   
            punchcardname=swoutdir+'sw.'+str(cpu3)+'.sh'   
            punchcard=open(punchcardname,'w+')
            punchcard.write(swbindir+'matrix_sw_vx <<EOF &> '
                            +logfil+'\n')
            punchcard.write('"'+swdatdir+'"\n')
            punchcard.write('"'+swoutdir+'"\n')
            punchcard.write('"'+keroutdir+'"\n')
            punchcard.write('"'+tmpnam+'"\n')
            punchcard.write(str(irorrl)+'\n')
            punchcard.write(matmode)
            punchcard.write(str(swdatun)+'\n')
            punchcard.write(str(eqincr)+'\n')
            punchcard.write(str(iadapt)+'\n')
            punchcard.write(str(nlay)+'\n')
            punchcard.write(str(nombloc)+'\n')
            punchcard.write(str(itera)+'\n')
            punchcard.write('"'+swdatset+'"\n')
            punchcard.write('EOF\n')
            punchcard.close()
            obs=obs.split()
            print 'Appending ',obs,' to submit script'
            nrofobs=line_count(swdatdir+'summary/summ.'+typ+'.'
                               +obs[1]+'.'+obs[2]+'.txt')
            os.system('chmod 755 '+punchcardname)       
            if int(nrofobs/ophr)<=7:
                que=7
            else:
                que=35
            if machine['sche']=='bsub':               
                submit.write('bsub -R lustre -W'+str(que)+':59 <'
                             +punchcardname+'\n')
            elif machine['sche']=='qsub':               
                submit.write('qsub -l nodes=1:ppn=1,walltime='+
                             str(que)+':00:00 '+punchcardname+'\n')
            elif machine['sche']=='pc':               
                submit.write(punchcardname+' &\n')

            cpu=cpu+1
            cpu3="%03d" % cpu
        submit.close()
        os.system('chmod 755 '+submitname)    


# ================================================================
#  Major subroutine that prepares submit scripts for body
#  wave matrix computations. This part of the code has not
#  been "parallelized" so far, but this would be trivial
#  and can be done if neeed
#  

def prepare_bowa_matrix(bwdatdir,bwoutdir,rootdir,eqincr,
                        bwdatset,layers,adaptive,xivoigt,
                        invpars):
    print 'Computing body wave matrices'
       
    if adaptive:
        iadapt=1
    else:
        iadapt=0
        
    bwbindir=rootdir+'/mat/bowa/' #default
    bwrefmod=bwbindir+'/refmods/aniprembw'

    # vph vpv vsh vsv part is kept general now
    matmode='"\n'+str(len(invpars))+'\n'
    bwbindir=rootdir+'/mat/bowa/'            
    for i in invpars:
        matmode=matmode+str(i)+'\n'
        
    # This part is still a separate fork
    if xivoigt:
        bwbindir=rootdir+'/mat/bowaav/'        
        matmode='"\n2\n3\n4\n'        
        
    obsfilbw=bwdatdir+'/observed_phases.txt'
    obsfilbw = obsfilbw.replace('//','/')
    bwdatdir = bwdatdir.replace('//','/')
    bwrefmod = bwrefmod.replace('//','/') 
    bwphasefile=open(obsfilbw)
    bwphaselist=bwphasefile.readlines()
    bwphaselist=[(el.strip('\n')) for el in bwphaselist]
    layfil = bwoutdir+'layers.in'
    layfil = layfil.replace('//','/')
    dump_layers(layers,layfil)
    for panda in bwphaselist:
        # this separates the phasename in form ScS.0
        # in arc and phase, 0=minor, 1=major        
        splitpanda=panda.split('.')
        phase=str(splitpanda[0])
        arc=str(splitpanda[1])
        print "Starting to generate matrices for "+phase
        print "Which has the arc-code: "+arc+" (0=minor; 1=major)"
        bwdatfil=bwdatdir+'database.'+panda
        punchcardname=bwoutdir+'bw.pc'
        punchcard=open(punchcardname,'w+')        
        # define command line string
        print bwoutdir
        writestr=[str(eqincr),'\n',
                  str(iadapt)+'\n"',
                  bwoutdir+'"\n"',
                  layfil,
                  matmode,
                  bwrefmod,'\n0\n1\n"',
                  phase,'"\n1\n0 360\n',
                  bwdatfil,'\n',str(arc)]        
        # pipe command into matrix_bw_vx
        for string in writestr:
            punchcard.write(string)
        punchcard.close()        
        # execute matrix_bw_vx
        cmdstr=bwbindir+'matrix_bw_vx < '+punchcardname
        os.system(cmdstr)
        # matrix files should have a nice name
        if arc=="1":
            phaseout=str(phase)+"m"
        elif arc=="0":
            phaseout=str(phase)
        else:
            print arc
            print "something went wrong"
            exit()
        # move matrix to correct folder
        appendixnew=['bw.xxx.'+bwdatset+'.', 'bw.indx.'+bwdatset+'.',
                     'bw.poin.'+bwdatset+'.', 'bw.rhs.'+bwdatset+'.']
        appendixold=['mat', 'ind', 'pnt', 'vec']
        for i in range(0,len(appendixnew)):
            src=bwoutdir+'/tmp/a.vx_'+appendixold[i]            
            des=bwoutdir+'/'+appendixnew[i]+phaseout+'.d'
            shutil.move(src,des)
        # move raypaths and sensitivity
        src=bwoutdir+'/tmp/sens.txt'
        des=bwoutdir+'sens_'+phaseout+'.txt'
        shutil.move(src,des)
        src=bwoutdir+'/tmp/ray.txt'
        des=bwoutdir+'ray_'+phaseout+'.txt'
        shutil.move(src,des)
        # move dws and htc
        src=bwoutdir+'/tmp/a.vx_dws'
        des=bwoutdir+'/dws/'+'dwsa.'+bwdatset+'.'+phaseout+'.d'
        shutil.move(src,des)

        # move actual source and receiver paths
        # needed for computation of corrections
        src=bwoutdir+'/tmp/quakes_actu.dat'
        des=bwoutdir+'quakes_actu_'+phaseout+'.dat'
        shutil.move(src,des)
        src=bwoutdir+'/tmp/receiv_actu.dat'
        des=bwoutdir+'receiv_actu_'+phaseout+'.dat'
        shutil.move(src,des)


# ================================================================
#  Major subroutine that prepares submit scripts for body
#  performing projections onto a data-adaptive grid.
#  

def prepare_projection(swfolders,bwfolders,adpoutdir,
                       phasearray,weigharray,swupscale,
                       bwupscale,synth,model,cfacsw,
                       cfacbw,rdamp,ddamp,ndamp,nlay,
                       eqincr,rootdir,thresh,olddws,
                       dwsfile,machine):

    # Some hardcoded parameters
    pcname = 'proj.'
    nonz      = int(9e8)     # Estimate of total nonzero entries
    regadfac  = int(5e5)     # Estimate of rows associated with regularization
    adpbindir = rootdir+'/adp/bin'
    submit    = open(adpoutdir+'/submit_2.sh','w+')
    iswit     = 1
    iadapt    = 1 # we are NOT in adaptive mode                  
    gridfin   = eqincr
    minhits   = 20
    outofeu   = 0

    # create upscale arrays
    availpaths=[]; availpathsdws=[]
    upscalearray=[]; upscalearraydws=[]    
    for item in swfolders:
        files = os.listdir(item.split()[0])
        for j in files:
            availpaths.append(item.split()[0]+'/'+j)
            upscalearray.append(swupscale[swfolders.index(item)])       
        files = os.listdir(item.split()[0]+'/dws/')
        for j in files:
            availpathsdws.append(item.split()[0]+'/dws/'+j)
            upscalearraydws.append(swupscale[swfolders.index(item)])
    for item in bwfolders:
        files = os.listdir(item.split()[0])
        for j in files:
            availpaths.append(item.split()[0]+'/'+j)
            upscalearray.append(bwupscale[bwfolders.index(item)])
        files = os.listdir(item.split()[0]+'/dws/')
        for j in files:
            availpathsdws.append(item.split()[0]+'/dws/'+j)
            upscalearraydws.append(bwupscale[bwfolders.index(item)])
                        

    # search absolute paths for patterns in phasearray
    count = 0
    templist = []
    for item in phasearray:
        for path in availpaths:
            if item in path:
                count=count+1
                templist.append(path)
                if count==4:
                    weigh1=float(weigharray[phasearray.index(item)])
                    weigh2=float(upscalearray[availpaths.index(path)])
                    floatweigh=weigh1*weigh2
                    weighting="%.3f" % floatweigh
                    templist.append(weighting)
                    count=0
    templist=sep_chunks(templist,5)

    # write list with phases for post-processing
    phaselist = []
    for phase in phasearray:
        for path in templist:
            for path2 in path:
                if phase in path2:
                    if '/sw.rhs' in path2 or '/bw.rhs.' in path2:
                        phaselist.append(phase)
    plist=open(adpoutdir+'/phaselist.log','w+')
    for i in phaselist:
        plist.write(str(i)+'\n')
    plist.close()
           

    # do the same search for the dws files, to get a
    # list of all actually used dws files and their
    # associated weight
    dwslist = []
    for item in phasearray:
        for path in availpathsdws:
            if item in path:
                dwslist.append(path)
                weigh1=float(weigharray[phasearray.index(item)])
                weigh2=float(upscalearraydws[availpathsdws.index(path)])
                floatweigh=weigh1*weigh2
                weighting= "%.3f" % floatweigh
                dwslist.append(weighting)           
    dwslist=sep_chunks(dwslist,2)

    # Consistency checks
    if len(dwslist)!=len(templist):
        print len(dwslist)
        print len(templist)
        print "*** ERROR ***"
        print "in prepare_project: It seems as if not enough .dws files are stored in"
        print "one of the matrix folders. Perhaps rerun the matrix submit script ... quitting!"
        print "*************"
        sys.exit()

    # another consistency check
    if len(templist)==0:
        print "*** ERROR ***"
        print "in prepare_project: Presumably you have correctly crawled one/some matrix folders"
        print "But weighted down all data to 0 using a wrong weighting scheme ... quitting!"
        print "*************"
        sys.exit()
    
    if size(templist)/len(templist)!=5:
        print "*** ERROR ***"
        print "in prepare_project: This error seems a bit more severe. Check whats going on!"
        print "*************"
        sys.exit()

    # Read dws files, add them up and write them out
    # to one file elone=dwslist[0]
    if not olddws:
        dws=np.zeros((line_count(dwslist[0][0]),4))
        dwsadd=open(adpoutdir+'/submit_1.sh','w+')
        dwsadd.write(adpbindir+'/addhtc <<EOF\n')        
        dwsadd.write(str(line_count(dwslist[0][0])/nlay)+'\n')
        dwsadd.write(str(nlay)+'\n')
        for item in dwslist:
            dwsadd.write('"'+item[0]+'"'+'\n')
            dwsadd.write(str(float(item[1]))+'\n')
        dwsadd.write('END\n')
        dwsadd.write('EOF\n')
        dwsadd.close()
        os.system('chmod 755 '+adpoutdir+'/submit_1.sh')        
    elif olddws:
        # copy existing dws file to working directory
        os.system('cp '+dwsfile+' '+adpoutdir+'/dws.all.dat')
    else:
        print 'Error in prepare_projection()'
        sys.exit()

    # put parts of matrix right order and write to pc, computes
    # cuttoff values and writes them to pc as well
    ntot=0; allident=[]
    matlist=open(adpoutdir+'/matrix.list.dat','w+')
    writegrid=0
    for item in templist:
        writegrid=writegrid+1
        # find the ident for the current file
        ident = dwslist[templist.index(item)][0].split('dws.')[1]
        allident.append(ident)
        # open punchcard for writing
        punchcard=open(adpoutdir+'/'+pcname+ident+'.sh','w+')
        punchcard.write(adpbindir+'/project <<EOF &> log.project.'+ident+'.dat\n')
        punchcard.write(str(iadapt)+'\n')    
        punchcard.write(str(nonz)+'\n')    
        punchcard.write(str(gridfin)+'\n')
        punchcard.write(str(iswit)+'\n')
        punchcard.write(str(outofeu)+'\n')
        punchcard.write(str(thresh[0])+'\n')    
        punchcard.write(str(thresh[1])+'\n')
        punchcard.write(str(minhits)+'\n')
        punchcard.write(str(nlay)+'\n')
        punchcard.write(str(int(writegrid))+'\n')
        bwind1=first_substr(item, '/bw.xxx.')
        bwind2=first_substr(item, '/bw.indx.')
        bwind3=first_substr(item, '/bw.poin.')
        bwind4=first_substr(item, '/bw.rhs.')
        swind1=first_substr(item, '/sw.xxx.')
        swind2=first_substr(item, '/sw.indx.')
        swind3=first_substr(item, '/sw.poin.')
        swind4=first_substr(item, '/sw.rhs.')  
        if bwind1 == -1 and bwind2 == -1 and bwind3 == -1 and bwind4 == -1:
            if swind1 != -1 and swind2 != -1 and swind3 != -1 and swind4 != -1:
                print 'writing out a sw to project pc'
                linenr=line_count(item[swind4])
                ntot=ntot+linenr
                punchcard.write('"'+item[swind1]+'"'+'\n')
                punchcard.write('"'+adpoutdir+'/sw.xxx.'+ident+'.d.ad"\n')                        
                punchcard.write('"'+item[swind2]+'"'+'\n')
                punchcard.write('"'+adpoutdir+'/sw.indx.'+ident+'.d.ad"\n')      
                punchcard.write('"'+item[swind3]+'"'+'\n')
                punchcard.write('"'+adpoutdir+'/sw.poin.'+ident+'.d.ad"\n')                        
                punchcard.write('"'+item[swind4]+'"'+'\n')
                punchcard.write('"'+adpoutdir+'/sw.rhs.'+ident+'.d.ad"\n')
                matlist.write('sw'+' '+ident+'.d.ad '+str(item[4])+'\n') 
            else:
                print 'Error in prepare_projection()'                
        if swind1 == -1 and swind2 == -1 and swind3 == -1 and swind4 == -1:
            if bwind1 != -1 and bwind2 != -1 and bwind3 != -1 and bwind4 != -1:
                print 'writing out a bw to project pc'
                linenr=line_count(item[bwind4])
                ntot=ntot+linenr
                punchcard.write('"'+item[bwind1]+'"'+'\n')
                punchcard.write('"'+adpoutdir+'/bw.xxx.'+ident+'.d.ad"\n')            
                punchcard.write('"'+item[bwind2]+'"'+'\n')
                punchcard.write('"'+adpoutdir+'/bw.indx.'+ident+'.d.ad"\n')                        
                punchcard.write('"'+item[bwind3]+'"'+'\n')
                punchcard.write('"'+adpoutdir+'/bw.poin.'+ident+'.d.ad"\n')                        
                punchcard.write('"'+item[bwind4]+'"'+'\n')
                punchcard.write('"'+adpoutdir+'/bw.rhs.'+ident+'.d.ad"\n')
                matlist.write('bw'+' '+ident+'.d.ad '+str(item[4])+'\n')
            else:
                print 'Error in prepare_projection()'                
        punchcard.write('"'+adpoutdir+'/dws.all.dat"\n')
        punchcard.write('"'+adpoutdir+'/grid.info.'+ident+'"\n')
        punchcard.write('"'+adpoutdir+'/adpx.info.'+ident+'"\n')
        punchcard.write('"'+adpoutdir+'/htca.info.'+ident+'"\n')
        punchcard.write('"'+adpoutdir+'/numa.info.'+ident+'"\n')
        punchcard.write(str(linenr)+'\n')
        punchcard.write('EOF\n')
        punchcard.close()
    matlist.close()
    print len(templist), ' subsets added to projection punchcard'    

    # reopen to write ntot and other parameters to each punchcard
    # and to generate the submit script
    submit=open(adpoutdir+'/submit_2.sh','w+')
    submit.write(adpoutdir+'/submit_1.sh\n')
    for item in allident:
        os.system('chmod 755 '+adpoutdir+'/'+pcname+item+'.sh')
        if machine["sche"]=='bsub':
            submit.write('bsub -R lustre -W'+str(machine['wall'])+':59 < '
                         +adpoutdir+pcname+item+'.sh\n')
        elif machine["sche"]=='qsub':
            submit.write('qsub -l nodes=1:ppn=1,walltime='+str(machine['wall'])+
                         ':00:00 '+adpoutdir+pcname+item+'.sh\n')
        elif machine["sche"]=='pc':
            submit.write(adpoutdir+pcname+item+'.sh &\n')

        punchcard.close()
    submit.close()

    # change submitscript to executable    
    os.system('chmod 755 '+adpoutdir+'submit_2.sh')   



# ================================================================
#  Major subroutine that prepares submit scripts for
#  regular grid inversions and many damping parameters
#  

def prepare_inversion_reg(swfolders,bwfolders,invoutdir,
                          phasearray,weigharray,swupscale,
                          bwupscale,synth,model,cfacsw,
                          cfacbw,rdamp,ddamp,ndamp,nlay,
                          eqincr,dmp_file_r,dmp_file_d,
                          rootdir,invpars,preani,preani_model,
                          machine):

    # Some parameters one may need to adapt
    nonz      = int(9e8)     # Estimate of total nonzero entries
    regadfac  = int(5e5)     # Estimate of rows associated with regularization
    invbindir = rootdir+'/inv/bin'
    sigma     = 3
    submit    = open(invoutdir+'/submit.sh','w+')
    iswit     = 1
    iadapt    = 0 # we are NOT in adaptive mode
    cpunum    = machine["ncpu"]
    time      = machine["wall"]
    nentries  = machine["nent"]
    memory    = machine["memo"]
    
    # create upscale arrays
    availpaths = []; upscalearray = []    
    for item in swfolders:
        files = os.listdir(item.split()[0])
        for j in files:
            availpaths.append(item.split()[0]+'/'+j)
            upscalearray.append(swupscale[swfolders.index(item)])       
    for item in bwfolders:
        files = os.listdir(item.split()[0])
        for j in files:
            availpaths.append(item.split()[0]+'/'+j)
            upscalearray.append(bwupscale[bwfolders.index(item)])

    # search paths for patterns in phasearray
    count = 0; templist = []
    for item in phasearray:
        for path in availpaths:
            if item in path:
                count=count+1
                templist.append(path)
                if count==4:
                    weigh1=float(weigharray[phasearray.index(item)])
                    weigh2=float(upscalearray[availpaths.index(path)])
                    floatweigh=weigh1*weigh2
                    weighting= "%.3f" % floatweigh
                    templist.append(weighting)
                    count=0
    templist=sep_chunks(templist, 5)

    # write list with phases for post-processing
    phaselist = []
    for phase in phasearray:
        for path in templist:
            for path2 in path:
                if phase in path2:
                    if '/sw.rhs' in path2 or '/bw.rhs.' in path2:
                        phaselist.append(phase)
    plist=open(invoutdir+'/phaselist.log','w+')
    for i in phaselist:
        plist.write(str(i)+'\n')
    plist.close()
       
    # another consistency check
    if len(templist)==0:
        print "*** ERROR ***"
        print "in prepare_inversion_reg: Presumably you have correctly crawled one/some matrix folders"
        print "But weighted down all data to 0 using a wrong weighting scheme ... quitting!"
        print "*************"
        sys.exit()
    
    if size(templist)/len(templist)!=5:
        print "*** ERROR ***"
        print "in prepare_inversion_reg: This error seems a bit more severe. Check whats going on!"
        print "*************"
        sys.exit()
        
    if synth:
        linemod=line_count(model)
        make_dir(invoutdir+'/synth/')
        synthpc=open(invoutdir+'/synthetics.sh','w+')
                
    # open punchcard for writing
    pcname0 = '/pc.tmp'
    punchcard=open(invoutdir+pcname0,'w+')    
    # put parts of matrix right order and write to pc, computes
    # cuttoff values and writes them to pc as well
    ntot=0

    for item in templist:
        bwi1=first_substr(item,'/bw.xxx.')
        bwi2=first_substr(item,'/bw.indx.')
        bwi3=first_substr(item,'/bw.poin.')
        bwi4=first_substr(item,'/bw.rhs.')        
        # synthetic or data inversion
        swi1=first_substr(item,'/sw.xxx.')
        swi2=first_substr(item,'/sw.indx.')
        swi3=first_substr(item,'/sw.poin.')
        swi4=first_substr(item, '/sw.rhs.')
        
        if bwi1==-1 and bwi2==-1 and bwi3==-1 and bwi4==-1:
            if swi1!=-1 and swi2!=-1 and swi3!=-1 and swi4!=-1:
                print 'Appending surface wave submatrix to punchcard'
                punchcard.write('"'+item[swi1]+'"'+'\n')
                punchcard.write('"'+item[swi2]+'"'+'\n')
                punchcard.write('"'+item[swi3]+'"'+'\n')
                print item[swi4]
                linenr=line_count(item[swi4])
                print str(linenr)+" rows"            
                if synth:
                    synthpc.write(invbindir+'/adotm <<EOF\n')
                    synthpc.write('"'+item[swi1]+'"'+'\n')
                    synthpc.write('"'+item[swi2]+'"'+'\n')
                    synthpc.write('"'+item[swi3]+'"'+'\n')
                    synthpc.write('"'+model+'"'+'\n')
                    synthpc.write(str(sigma)+'\n')
                    idsynrhs=item[swi4].split('/sw.rhs.')[1]
                    synthpc.write('"'+invoutdir+'/synth/sw.rhs.'+idsynrhs+'"\n')
                    synthpc.write(str(linemod)+'\n')
                    synthpc.write(str(linenr)+'\n')
                    synthpc.write('EOF\n')
                    punchcard.write('"'+invoutdir+'/synth/sw.rhs.'+idsynrhs+'"\n')
                else:
                    punchcard.write('"'+item[swi4]+'"'+'\n')                                   
                ntot=ntot+linenr            
                # Compute cutoff and conv to str
                data=np.genfromtxt(item[swi4], dtype=None)
                cutoffu=np.mean(data)+(cfacsw*np.std(data))
                cutoffl=np.mean(data)-(cfacsw*np.std(data))
                cutoffu="%.3f" % cutoffu
                cutoffl="%.3f" % cutoffl            
                # Write extra info for dataset
                punchcard.write(item[4]+'\n')       
                punchcard.write(cutoffu+'\n')
                punchcard.write(cutoffl+'\n')
                punchcard.write(str(linenr)+'\n')
            else:
                print 'Error in prepare_matrix_reg()'
                sys.exti()

        if swi1==-1 and swi2==-1 and swi3==-1 and swi4==-1:
            if bwi1!=-1 and bwi2!=-1 and bwi3!=-1 and bwi4!=-1:
                print 'Appending body wave matrix to punchcard'                
                punchcard.write('"'+item[bwi1]+'"'+'\n')
                punchcard.write('"'+item[bwi2]+'"'+'\n')
                punchcard.write('"'+item[bwi3]+'"'+'\n')
                print item[bwi4]
                linenr=line_count(item[bwi4])                
                print str(linenr)+" rows"
                if synth:
                    synthpc.write(invbindir+'/adotm <<EOF\n')
                    synthpc.write('"'+item[bwi1]+'"'+'\n')
                    synthpc.write('"'+item[bwi2]+'"'+'\n')
                    synthpc.write('"'+item[bwi3]+'"'+'\n')
                    synthpc.write('"'+model+'"'+'\n')
                    synthpc.write(str(sigma)+'\n')
                    idsynrhs=item[bwi4].split('/bw.rhs.')[1]
                    synthpc.write('"'+invoutdir+'/synth/bw.rhs.'+idsynrhs+'"\n')
                    synthpc.write(str(linemod)+'\n')
                    synthpc.write(str(linenr)+'\n')
                    synthpc.write('EOF\n')   
                    punchcard.write('"'+invoutdir+'/synth/bw.rhs.'+idsynrhs+'"\n')                    
                else:
                    punchcard.write('"'+item[bwi4]+'"'+'\n')
                ntot=ntot+linenr
                data=np.genfromtxt(item[bwi4], dtype=None)
                # Compute cutoff and conv to str
                cutoffu=np.mean(data)+(cfacbw*np.std(data))
                cutoffl=np.mean(data)-(cfacbw*np.std(data))
                cutoffu="%.3f" % cutoffu
                cutoffl="%.3f" % cutoffl            
                # Write extra info for dataset
                punchcard.write(item[4]+'\n')
                punchcard.write(cutoffu+'\n')
                punchcard.write(cutoffl+'\n')
                punchcard.write(str(linenr)+'\n')
            else:
                print 'Error in prepare_matrix_reg()'
                sys.exit()

    if synth:        
        synthpc.close()
        os.system('chmod 755 '+invoutdir+'/synthetics.sh')
        os.system(invoutdir+'/synthetics.sh')
        print "DONE WITH SYNTHETICS!"
           

    print len(templist), ' subsets added to punchcard'
    punchcard.write('finished\n'); punchcard.close()    
    punchcard0=open(invoutdir+pcname0,'r+')
    oldpc0 = punchcard0.read() # read everything in the file
    ntot=ntot+regadfac 
    nonz=ntot*nentries
    kk=0

    # Loop over all inversion runs
    for i in range(0,len(rdamp)):
        kk=kk+1
        run3="%03d" % kk
        runid=str(rdamp[i])+'.'+str(ddamp[i])+'.'+\
              str(int(nlay))+'.'+str(int(eqincr))
        pcname = 'pc.'+runid+'.sh'
        punchcard=open(invoutdir+pcname,'w+')
        punchcard.write(oldpc0)   # write the new line before
        punchcard.write('"'+dmp_file_r+'"\n')
        punchcard.write('"'+dmp_file_d+'"\n')
        punchcard.write(str(rdamp[i])+'\n')       
        punchcard.write(str(ddamp[i])+'\n')
        punchcard.write(str(ndamp[i])+'\n')        
        if preani>0.0 or preani==-1.0:
            punchcard.write(str(preani)+'\n')
            punchcard.write('"'+preani_model+'"\n')
        else:
            punchcard.write(str(preani)+'\n')                        
        punchcard.write('"'+runid+'"\n')        
        punchcard.close()
        # Insert matrix bytelengths, etc
        with open(invoutdir+pcname, "r+") as punchcard:
            oldpc = punchcard.read()   # re-read the file
            punchcard.seek(0) # rewind
            punchcard.write(invbindir+'/flexinv90 <<EOF &>log.'+runid+'.dat\n')
            punchcard.write(str(ntot)+'\n')
            punchcard.write(str(nonz)+'\n')
            punchcard.write(str(len(invpars))+'\n')
            punchcard.write(str(eqincr)+'\n')
            punchcard.write(str(iswit)+'\n')
            punchcard.write(str(nlay)+'\n')
            punchcard.write('"'+str(invoutdir)+'"'+'\n')
            punchcard.write(str(iadapt)+'\n')
            punchcard.write(oldpc)     # write the new line before
            punchcard.write('EOF\n')   # write the new line before
        punchcard.close()
        if machine["sche"]=='bsub':
            submit.write('bsub -n '+str(cpunum)+' -R "rusage[mem='+str(memory)+\
                         ']" -R lustre -W'+str(time)+':59 < '\
                         +invoutdir+pcname+'\n')
        elif machine["sche"]=='qsub':
            submit.write('qsub -l nodes=1:ppn=1,mem='+str(memory)+'mb,walltime='\
                         +str(time)+':00:00 '+invoutdir+pcname+'\n')
        elif machine["sche"]=='pc':
            submit.write(invoutdir+pcname+' &\nwait\n')

    submit.close()
    os.system('chmod 755 '+\
              invoutdir+'/submit.sh')   


# ================================================================
#  Major subroutine that prepares submit scripts for
#  adaptive grid inversions and many damping parameters
#  

def prepare_inversion_adp(adpindir,invoutdir,synth,model,cfacsw,cfacbw,rdamp,
                          ddamp,ndamp,nlay,eqincr,dmp_file_r,dmp_file_d,rootdir,
                          invpars,machine):

    # Some hardcoded parameters
    nonz      = int(9e8)     # Estimate of total nonzero entries
    regadfac  = int(5e5)     # Estimate of rows associated with regularization
    invbindir = rootdir+'/inv/bin'
    submit    = open(invoutdir+'/submit.sh','w+')
    sigma     = 3
    iswit     = 1
    iadapt    = 1 # we are in adaptive mode
    time      = machine["wall"]
    nentries  = machine["nent"]
    memory    = machine["memo"]
    cpunum    = machine["ncpu"]
    
    print "Performing an inversion with an adaptive grid"
    print "Combining all data which is contained in the folder ",adpindir[0]

    # Read the list which was created by prepare_project.py
    matadap=[]
    matlist=open(adpindir[0]+'/matrix.list.dat')
    matadap=matlist.readlines()
    
    # check how many adaptive gridpoints
    grdinf=glob.glob(adpindir[0]+'/grid.info.*.d')[0]
    npixadapt=line_count(grdinf)

    # synthetics?
    if synth:
        linemod=line_count(model)
        make_dir(invoutdir+'/synth/')
        synthpc=open(invoutdir+'/synthetics.sh','w+')

    # open punchcard for writing
    pcname0 = '/pc.sh'
    punchcard=open(invoutdir+pcname0,'w+')
    # put parts of matrix right order and write to pc, computes
    # cuttoff values and writes them to pc as well
    ntot=0; templist=[]
    for item0 in matadap:
        item=item0.split()
        templist.append("."+item[1]+".")
        punchcard.write('"'+adpindir[0]+'/'+item[0]+'.xxx.'+item[1]+'"\n')
        punchcard.write('"'+adpindir[0]+'/'+item[0]+'.indx.'+item[1]+'"\n')
        punchcard.write('"'+adpindir[0]+'/'+item[0]+'.poin.'+item[1]+'"\n')
        rhsname=adpindir[0]+'/'+item[0]+'.rhs.'+item[1]
        linenr=line_count(rhsname)
        if synth:
            synthpc.write(invbindir+'/adotm <<EOF\n')
            synthpc.write('"'+adpindir[0]+'/'+item[0]+'.xxx.'+item[1]+'"\n')
            synthpc.write('"'+adpindir[0]+'/'+item[0]+'.indx.'+item[1]+'"\n')
            synthpc.write('"'+adpindir[0]+'/'+item[0]+'.poin.'+item[1]+'"\n')
            synthpc.write('"'+model+'"'+'\n')
            synthpc.write(str(sigma)+'\n')
            idsynrhs=item[1] #.split('/sw.rhs.')[1]
            synthpc.write('"'+invoutdir+'/synth/sw.rhs.'+idsynrhs+'"\n')
            synthpc.write(str(linemod)+'\n')
            synthpc.write(str(linenr)+'\n')
            synthpc.write('EOF\n')
            punchcard.write('"'+invoutdir+'/synth/sw.rhs.'+idsynrhs+'"\n')
        else:
            punchcard.write('"'+adpindir[0]+'/'+item[0]+'.rhs.'+item[1]+'"\n')
        ntot=ntot+linenr            
        # Compute cutoff and conv to str
        data=np.genfromtxt(rhsname, dtype=None)
        if item[0]=='sw':
            print 'writing surface wave matrix path to punchcard'
            cfac=cfacsw
        elif item[0]=='bw':
            print 'writing body wave matrix path to punchcard'
            cfac=cfacbw
        cutoffu=np.mean(data)+(cfac*np.std(data))
        cutoffl=np.mean(data)-(cfac*np.std(data))
        cutoffu="%.3f" % cutoffu
        cutoffl="%.3f" % cutoffl            
        # Write extra info for dataset
        punchcard.write(item[2]+'\n')       
        punchcard.write(cutoffu+'\n')
        punchcard.write(cutoffl+'\n')
        punchcard.write(str(linenr)+'\n')
    print len(matadap), ' datasets appended to punchcard'
    punchcard.write('finished\n')
    punchcard.close()
    punchcard0=open(invoutdir+pcname0,'r+')
    oldpc0 = punchcard0.read() # read everything in the file
    ntot=ntot+regadfac 
    nonz=ntot*nentries
    kk=0

    if synth:        
        synthpc.close()
        os.system('chmod 755 '+invoutdir+'/synthetics.sh')
        os.system(invoutdir+'/synthetics.sh')
        print "DONE WITH SYNTHETICS!"

    # Loop over all inversion runs
    for i in range(0,len(rdamp)):
        kk=kk+1
        run3="%03d" % kk
        runid=str(rdamp[i])+'.'+str(ndamp[i])+'.'+\
              str(int(nlay))+'.'+str(int(eqincr))        
        pcname = 'pc.'+runid+'.sh'
        punchcard=open(invoutdir+pcname,'w+')
        punchcard.write(oldpc0)   # write the new line before
        punchcard.write('"'+dmp_file_r+'"\n')
        punchcard.write('"'+dmp_file_d+'"\n')        
        punchcard.write(str(rdamp[i])+'\n')
        punchcard.write(str(ddamp[i])+'\n')
        punchcard.write(str(ndamp[i])+'\n')        
        punchcard.write('"'+runid+'"\n')        
        punchcard.close()
        # Insert matrix bytelengths, etc
        with open(invoutdir+pcname, "r+") as punchcard:
            oldpc = punchcard.read() # read everything in the file
            punchcard.seek(0) # rewinds
            punchcard.write(invbindir+'/flexinv90 <<EOF &> log.'\
                            +runid+'.dat\n')
            punchcard.write(str(ntot)+'\n')
            punchcard.write(str(nonz)+'\n')
            punchcard.write(str(len(invpars))+'\n')
            punchcard.write(str(eqincr)+'\n')
            punchcard.write(str(iswit)+'\n')
            punchcard.write(str(nlay)+'\n')
            punchcard.write('"'+str(invoutdir)+'"'+'\n')
            punchcard.write(str(iadapt)+'\n')
            punchcard.write('"'+str(grdinf)+'"'+'\n')
            punchcard.write(oldpc)   # write the new line before
            punchcard.write('EOF\n')   # write the new line before
        punchcard.close()
        if machine["sche"]=='bsub':
            submit.write('bsub -n '+str(cpunum)+' -R "rusage[mem='+str(memory)+\
                         ']" -R lustre -W'+str(time)+':59 < '\
                         +invoutdir+pcname+'\n')
        elif machine["sche"]=='qsub':
            submit.write('qsub -l nodes=1:ppn=1,mem='+str(memory)+'mb,walltime='\
                         +str(time)+':00:00 '+invoutdir+pcname+'\n')
        elif machine["sche"]=='pc':
            submit.write(invoutdir+pcname+' &\nwait\n')
    submit.close()
    os.system('cp '+adpindir[0]+'/phaselist.log '+invoutdir)
    os.system('cp '+str(grdinf)+'.lay '+invoutdir+'/grid.info.lay')
    os.system('cp '+str(grdinf)+'.sh '+invoutdir+'/grid.info.sh')
    os.system('chmod 755 '+invoutdir+'/submit.sh')   


# ================================================================
#  call subroutine main
#  

if __name__ == "__main__":

    main(sys.argv[1:])
