# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:50:07 2018

@author: M.C. Jose Luis Rodriguez-Solis
         Departamento de Oceanografia Fisica
         Centro de Investigación Científica y de Educación Superior de Ensenada
         CICESE

Este programa es para hacer promedios mensuales con alguna variable
de los datos del modelo CESM-AMIP

"""

from netCDF4 import Dataset, date2num, num2date

import numpy as np
import scipy.signal as signal

import time as pytime
from calendar import monthrange, isleap
from dateutil.relativedelta import relativedelta
import datetime
tic = pytime.clock()
import os
import getpass
import glob
LinuxUser = getpass.getuser()
# Caso
exp_case='camip'
varname='T'
varnum ='01'

listadata=glob.glob('/scratch/201803026n-2/cesm/var/'+exp_case+'/var'+varnum+'/*.nc')
listadata.sort()

dt = (datetime.date(2005, 12, 31) - datetime.date(1981, 1, 1)).days
# leap years count
yleap = 0
for i in range(1981,2006):
    if isleap(i):
        yleap = yleap + 1
dt = dt - yleap + 15
# extrayendo algunas variables
nc       = Dataset (listadata[0])
timed    = nc.variables['time'][:]
unitsd   = nc.variables['time'].units
x        = nc.variables['lon'][:]
y        = nc.variables['lat'][:]
nlev     = nc.variables['lev'][:]
longname = nc.variables[varname].long_name
VARunits = nc.variables[varname].units

[nt, nz, ny, nx] = nc.variables[varname].shape
[xd,yd] = np.meshgrid (x,y)

year = datetime.datetime(1981,1,1).year

varmean = np.zeros((12*25,nz,ny,nx))

si = datetime.datetime(1980,12,31)
sf = datetime.datetime(2006,1,1)
for k in range (0,nz):
    nm = 0
    #uc = np.zeros((dt,ny,nx))* np.nan
    varc = np.zeros((dt,ny,nx))* np.nan
    ni = 0
    for n in range(0,len(listadata)):
        nc = Dataset (listadata[n])
        timed = nc.variables['time'][:]
        syear = num2date(timed,units=unitsd)
        fi = np.where(syear>si)
        nc.close() 
        try:
            if syear[fi[0][0]] > si and syear[fi[0][0]] < sf:
                 fi = fi[0][0]
                 print (syear[fi])
                 nc = Dataset (listadata[n])
                 var = nc.variables[varname][fi::,k,:,:].squeeze()
                 [nt, ny, nx] = var.shape
                 varc[ni:ni+nt,:,:] = var[:]
                 nc.close()
                 del  var
                 print (ni,nt)
                 ni = ni + nt
        except:
            print (syear[0])
            print ('no entro')

    [nt,ny,nx] = varc.shape
    print (n)
    a = datetime.date(year,1,1)
    idi = 0
    ti = []
    for t in range (0,nt,1):
        lm = monthrange(a.year,a.month)[1]
        if lm == 29:
            lm = 28
        idf = idi + lm
        # promedio mensual
        varvt = np.nanmean(varc[idi:idf,:,:], axis=0)
        varmean[nm,k,:,:] = varvt

        print ('>> month number :' + str(nm) + ' ' + a.strftime('%Y-%m'))

        nm = nm+1

        if k == 29:
            b  = datetime.datetime(a.year,a.month,a.day)
            b = date2num(b,unitsd,calendar='standard')
            ti.append(b)
        a = a + relativedelta(months=1)
        print (t,idi,idf,nt)
        idi = idf*1
        if a.year == 2006 and a.month == 01:
            break
         
    print (a.strftime('%Y-%m'))
    print ()
    print ('level ' + str(nlev[k]) + ' mb Done... ')
    print ()
ti = np.asarray(ti)

print ('-----> Saving data in netcdf file...')
print ('')
ncname = varname+'_m_mean_'+exp_case+'.nc'
os.system("rm {0}".format(ncname))
# abriendo el archivo netcdf
f = Dataset(ncname,'w', format='NETCDF4') #'w' stands for write

## creando las dimenciones
f.createDimension('lon', len(x))
f.createDimension('lat', len(y))
f.createDimension('lev', len(nlev))
f.createDimension('time', None)

## creando las variables 
lon     = f.createVariable('lon', 'f8', 'lon')                         #double
lat     = f.createVariable('lat', 'f8', 'lat')                         #double
lev     = f.createVariable('lev', 'f8', 'lev')                         #double
time    = f.createVariable('time', 'f8', 'time', fill_value = False)   #double

var = f.createVariable(varname, 'f4', ('time','lev','lat', 'lon')) #float

## long names
lon.long_name = "longitude"
lat.long_name = "latitude"
lev.long_name = "pressure_level"

var.long_name = longname

## units
lon.units  = "degrees_east"
lat.units  = "degrees_north"
lev.units  = "millibars"
time.units = unitsd
var.units  = VARunits

## time calendar
time.calendar = "standard"

## asignando valores
lon[:]  = x[:]
lat[:]  = y[:]
lev[:]  = nlev[:]
time[:] = ti[:]
var[:]  = varmean[:]

f.history = "calculos con datos de experimento "+exp_case

f.close()

print ()
print ('Program done... you can check now the mat file with Eddy Momentum Flux Convergence [EMFC]')
toc = pytime.clock()
print
print ('time used in this process '+ str((toc -tic)/(60*60)) + ' hrs')
print ()
print ('Bye... ')
                         
