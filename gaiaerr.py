#!/usr/bin/env python                    
import numpy as np
import scipy
import sys
import os
import os.path
import argparse
import myutils

#--------Version History----------------------------------------------------------------------------
# 16/nov/2017: parallax mission time scaling factor fixed (goes as (5./tm)**0.5 not **1.)
# 11/oct/2016: VX,VY,VZ unit error fixed (inputs must be passes in mas/yr always, not muas/yr)

#Gaia error code path
gerr_path='/workd/cmateu/gaia_errors_color_tmission'
#gerr_path='/Users/cmateu/trabajo/gaia/gaia_challenge2014_mgc3/gaiaerror_py/'+'gaia_errors_color_tmission'

parser = argparse.ArgumentParser(description='Simulate Gaia errors + constant relative error in distance')
parser.add_argument('infile',metavar='infile(.ne.dat)',help='Input File (x y z vx vy vz Mv VI)',nargs=1,action='store')
parser.add_argument('-tm','--mission_t',help='Gaia mission time span in yr. Default 5.', action='store',default=5.,type=np.float)

parser.add_argument('-v','--verbose',help='Verbose', action='store_true',default=False)

#parse arguments
args = parser.parse_args()
infilen=args.infile[0]
mission_t=args.mission_t

if args.verbose: 
  print 'Input file:', infilen
  print 'Gaia Mission time:',mission_t

#Compute error scaling factor based on mission time (following Brown and deBruijne's prescriptions, priv. comm.)
if mission_t<=10.:
 pfactor=(5./mission_t)**0.5  #nominal errors are for mission_t=5., so factor==1 in this case. For parallax err scales as t
 factor=(5./mission_t)**1.5   #nominal errors are for mission_t=5., so factor==1 in this case
else:
 pfactor=(5./mission_t)**0.5   #If new Gaia is launched, scaling can be conservatively assumed to go as t
 factor=(5./mission_t)**1.    #If new Gaia is launched, scaling can be conservatively assumed to go as t
#Extra labels
if mission_t==5: tlabel=''
else: tlabel='%.1f' % (mission_t)

#Print auxiliary input file for gaerr code.Aux files have to be unique so multiple threads can be run simultaneuosly
auxinf=infilen+tlabel+'.aux.in'
auxoutf=infilen+tlabel+'.aux.out'
auxfilein=open(auxinf,'w')
auxfilein.write('%s\n%s\n' % (infilen,auxoutf))
auxfilein.close()

#Check whether auxfiles used by gaiaerr code exist in the present dir. Create symbolic links if not
if not os.path.isfile('avdisk.dat'):
  if args.verbose: print 'Gaia error code aux files missing, creating symbolic links...'
  proc='ln -s %s/*.dat .' % (gerr_path)
  os.system(proc)

#Run Gaia error code
if args.verbose: print 'Running Gaia error code...' 
os.system('%s/compute_err_color_gaia_tmission < %s' % (gerr_path,auxinf))

#Read gaia error output file
dat=scipy.genfromtxt(auxoutf)

#Get true parallax, simulate gpar by adding gaussian X% error 
relerr_par=dat[:,5-1]
xpar=dat[:,12-1]
gpar=dat[:,25-1]

#Rescale parallax err
sigma_par_new=(gpar-xpar)*pfactor
gpar=xpar+sigma_par_new
sigma_par_ref=relerr_par*xpar  #sigma used to draw random gaussian par error
sigma_par_ref_new=sigma_par_ref*pfactor
relerr_par_obs_new=sigma_par_ref_new/gpar   #this is fobs. relerr_par is ftrue
relerr_par_tru_new=sigma_par_ref_new/xpar

#Recompute gvrad (a lot come out from Merce's code as ****)
xvrad=dat[:,18-1]
sigma_vrad=dat[:,34-1]
gvrad=xvrad+np.random.normal(loc=0.,scale=sigma_vrad,size=xvrad.size)

gl,gb,gmulstar,gmub=dat[:,26-1],dat[:,27-1],dat[:,29-1],dat[:,30-1]
xl,xb,xmulstar,xmub=dat[:,13-1],dat[:,14-1],dat[:,16-1],dat[:,17-1]
#Recompute uncertainties
sigma_mulstar_new=(gmulstar-xmulstar)*factor
sigma_mub_new=(gmub-xmub)*factor
#Recompute 'observed proper motions'
gmulstar=xmulstar+sigma_mulstar_new
gmub=xmub+sigma_mub_new

fp=1000.
#Parallax for my function must be in muas. Proper motions must be in mas/yr (as needed by bovy library)
mydat=myutils.helio_obj(gl,gb,fp*gpar,gmulstar,gmub,gvrad,degree=True,flag_mulstar=True)

#Replace cols appropiately in full matrix
dat[:,25-1]=gpar
dat[:,26-1]=gl
dat[:,27-1]=gb
dat[:,28-1]=mydat.Rhel
dat[:,29-1]=gmulstar
dat[:,30-1]=gmub
dat[:,31-1]=gvrad
#---Proper motion cols----
dat[:,33-1]=dat[:,33-1]*factor  #sigma_mub
dat[:,36-1]=dat[:,36-1]*factor  #relerr_mub
#Parallax error cols------
dat[:, 5-1]=relerr_par_tru_new  #relerr_par
dat[:,32-1]=relerr_par_obs_new  #relerr_par_obs
#---Cartesian coords
dat[:,19-1]=-mydat.x
dat[:,20-1]=mydat.y
dat[:,21-1]=mydat.z
dat[:,22-1]=-mydat.vx
dat[:,23-1]=mydat.vy
dat[:,24-1]=mydat.vz

#Header and print formats
head_l=['Av','xV','Gmag','Grvs','relerr_par','xX','xY','xZ','xVX','xVY','xVZ','xpar_mas','xl_deg','xb_deg','xRhel','xmuls_cosb_mas','xmub_mas','xvrad','gX','gY','gZ','gVX','gVY','gVZ','gpar_mas','gl_deg','gb_deg','gRhel','gmuls_cosb_mas','gmub_mas','gvrad','relerr_parobs','sig_mub','sig_vrad','VI','relerr_mub','relerr_vrad']
head_cols=np.arange(len(head_l))+1
hfmts='#%17s '+(len(head_l)-1)*'%18s '
hfmts=hfmts+'\n'
fmts=(dat[0,:].size)*'%18.10f '

#Final output file name
ofilen=infilen.replace('.ne.dat','')+'.ge'+tlabel+'.dat'

#Print output file
if args.verbose: print 'Printing outputfile',ofilen
ofile=open(ofilen,'w')
ofile.write('#Gaia mission time assumed %.1f yr, error scaling factor %.3f\n' % (mission_t,factor))
ofile.write(hfmts % tuple(head_cols))
ofile.write(hfmts % tuple(head_l))
scipy.savetxt(ofile,dat,fmt=fmts)

#Remove aux files
proc='rm -f %s %s' % (auxinf,auxoutf)
os.system(proc)

#proc='rm -f TableVr-Jun2015.dat avdloc.dat avspir.dat myfile.ne.dat allruns.dat avori.dat gfactor-Jun2013.dat rf_allsky.dat avdisk.dat avori2.dat myfile.ge.dat run.dat'

if args.verbose: print 'Done'

