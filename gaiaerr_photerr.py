#!/usr/bin/env python                    
import numpy as np
import scipy
import sys
import os
import argparse
import myutils

#Gaia error code path
gerr_path='/Users/cmateu/trabajo/gaia/gaia_challenge2014_mgc3/gaiaerror_py/'+'gaia_errors_color_tmission'


parser = argparse.ArgumentParser(description='Simulate Gaia errors + constant relative error in distance')
parser.add_argument('infile',metavar='infile(.ne.dat)',help='Input File (x y z vx vy vz Mv VI)',nargs=1,action='store')
parser.add_argument('relerr_par',metavar='relative_par_error',help='Relative parallax error (constant)',nargs=1,action='store',type=np.float)
parser.add_argument('-tm','--mission_t',help='Gaia mission time span in yr. Default 5.', action='store',default=5.,type=np.float)

parser.add_argument('-v','--verbose',help='Verbose', action='store_true',default=False)

#parse arguments
args = parser.parse_args()
infilen=args.infile[0]
relerr_par=args.relerr_par[0]
mission_t=args.mission_t
if relerr_par>1:
  sys.exit('Relative Parallax Error larger than 100%... exiting.')
if args.verbose: 
  print 'Input file:', infilen
  print 'Relative Parallax error:',relerr_par
  print 'Gaia Mission time:',mission_t

#Compute error scaling factor based on mission time (following Brown and deBruijne's prescriptions, priv. comm.)
if mission_t<=10.:
 factor=(5./mission_t)**1.5   #nominal errors are for mission_t=5., so factor==1 in this case
else:
 factor=(5./mission_t)**1.    #If new Gaia is launched, scaling can be conservatively assumed to go as t
#Extra labels
if mission_t==5: tlabel=''
else: tlabel='%.0f' % (mission_t)

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
xpar=dat[:,12-1]
#gaiapar=gpar
gpar=xpar + xpar*np.random.normal(loc=0.,scale=relerr_par,size=xpar.size)

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
#Recompute mub_relerr
#relerr_mub=np.abs(sigma_mub_new/xmub)
#Recompute observed l,b
gl=xl+(gl-xl)*factor
gb=xb+(gb-xb)*factor

fp,fm=1000.,1000.
#Inputs for my function must be in muas
mydat=myutils.helio_obj(gl,gb,fp*gpar,fm*gmulstar,fm*gmub,gvrad,degree=True,flag_mulstar=True)

#Replace cols appropiately in full matrix
dat[:,25-1]=gpar
dat[:,26-1]=gl
dat[:,27-1]=gb
dat[:,28-1]=mydat.Rhel
dat[:,29-1]=gmulstar
dat[:,30-1]=gmub
dat[:,31-1]=gvrad
#---Rel err cols----
dat[:, 5-1]=relerr_par
dat[:,32-1]=(mydat.Rhel-dat[:,15-1])/dat[:,15-1] #(gRhel-xRhel)/xRhel
dat[:,33-1]=dat[:,33-1]*factor  #sigma_mub
dat[:,36-1]=dat[:,36-1]*factor  #relerr_mub
#---Cartesian coords
dat[:,19-1]=-mydat.x  #- so it matches the transformation used in Merce's code
dat[:,20-1]=mydat.y
dat[:,21-1]=mydat.z
dat[:,22-1]=mydat.vx
dat[:,23-1]=mydat.vy
dat[:,24-1]=mydat.vz

#Header and print formats
head_l=['Av','xV','Gmag','Grvs','relerr_par','xX','xY','xZ','xVX','xVY','xVZ','xpar_mas','xl_deg','xb_deg','xRhel','xmuls_cosb_mas','xmub_mas','xvrad','gX','gY','gZ','gVX','gVY','gVZ','gpar_mas','gl_deg','gb_deg','gRhel','gmuls_cosb_mas','gmub_mas','gvrad','relerr_D','sig_mub','sig_vrad','VI','relerr_mub','relerr_vrad']
head_cols=np.arange(len(head_l))+1
hfmts='#%17s '+(len(head_l)-1)*'%18s '
hfmts=hfmts+'\n'
fmts=(dat[0,:].size)*'%18.10f '

#Final output file name
ofilen=infilen.replace('.ne.dat','')+'.pe'+tlabel+'.dat'

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

#Gzip output file
#proc='gzip -f %s' % (ofilen)
#os.system(proc)

if args.verbose: print 'Done'

