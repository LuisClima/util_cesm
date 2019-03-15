# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:50:07 2018

@author: M.C. Jose Luis Rodriguez-Solis
         Departamento de Oceanografia Fisica
         Centro de Investigación Científica y de Educación Superior de Ensenada
         CICESE

Este programa lee todos los archivos del ERA Interim obtenidos del modelo CESM, con ellos
se calcula primero un filtro pasa altas para cada nivel y despues  se calcula  la
convergencia de flujo de momento 

Eddy Momemtum Flux Convergence [EMFC]


"""

from netCDF4 import Dataset, date2num, num2date

import numpy as np
import scipy.signal as signal
#import scipy.io as sio

import time as pytime
from calendar import monthrange, isleap
from dateutil.relativedelta import relativedelta
import datetime
tic = pytime.clock()
#import sys
import os
import getpass
import glob
LinuxUser = getpass.getuser()
# Caso
exp_case='lamip'

listadata01=glob.glob('/scratch/201803026n-2/cesm/var/'+exp_case+'/var00/*.nc')
listadata01.sort()
listadata02=glob.glob('/scratch/201803026n-2/cesm/var/'+exp_case+'/var01/*.nc')
listadata02.sort()

# Datos para el filtro
# la ventana en dias
ndays = 6
# el grado del polinomio
N = 5
# frecuencia muestral
fs = 1/(ndays*1.0)    # dt/myband in days
nl = 13
# filtro
filtro = 'highpass'
# spin up, n in years
ny = 2
spinup = (ny * 12)
#spinup = 2
# coeficientes del polinomio
B, A = signal.butter(N, fs, btype=filtro, output='ba')
# datos para la derivada meridional
R    = 6371000.0                 #Radio de la tierra
rad  = np.pi/180                 #lat a radianes
# length period
dt = (datetime.date(2005, 12, 31) - datetime.date(1981, 1, 1)).days
# leap years count
yleap = 0
for i in range(1981,2006):
    if isleap(i):
        yleap = yleap + 1
dt = dt - yleap + 1
# extrayendo algunas variables
nc = Dataset (listadata01[0])
names = tuple(nc.variables.keys())
timed     = nc.variables['time'][:]
unitsd    = nc.variables['time'].units
#calendard    = nc.variables['time'].calendar
#datesec  = nc.variables['datesec'][:]
#dated    = nc.variables['date'][:]
x        = nc.variables['lon'][:]
y        = nc.variables['lat'][:]
nlev     = nc.variables['lev'][:]
[nt, nz, ny, nx] = nc.variables['U'].shape
[xd,yd] = np.meshgrid (x,y)
year = datetime.datetime(1981,1,1).year
# datos para la derivada meridional
cos2 = (np.cos(yd * rad))**2
cos1 = np.cos(yd * rad)

print ("There are: "+ str(len(listadata01)) + " monthly CAM files")
print ()
print (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
print ("CAREFUL!!!")
print ("-----> The initial year file should be : " + listadata01[0][-39::])
print ("If file is correct, continue. There will be a 0 seconds pause")
print ("\33[32mPAUSE\033[0m'")
print (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
pytime.sleep(0)
print ()

FluxWind = np.zeros((12*25,nz,ny,nx))
MomeWind = np.zeros((12*25,nz,ny,nx))
Fluxtemp = np.zeros((12*25,nz,ny,nx))

print ()

si = datetime.datetime(1980,12,31)

for k in range (0,nz):
    nm = 0
    uc = np.zeros((dt,ny,nx))* np.nan
    vc = np.zeros((dt,ny,nx))* np.nan
    tc = np.zeros((dt,ny,nx))* np.nan
    print (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    print ('Begins level ' + str(nlev[k]) + ' mb, corresponding to: ' +str(k+1)+ ' of '+str(nz) )
    print (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    print ()
    ni = 0
    for n in range(0,len(listadata01)):
        #print ('reading > '+listadata01[n][-39::],n+1)
        nc = Dataset (listadata01[n])
        timed = nc.variables['time'][:]
        syear = num2date(timed,units=unitsd)
        fi = np.where(syear>si)
        try:
            if syear[fi[0][0]] > si:
                 fi = fi[0][0]
                 print (syear[fi])
                 u = nc.variables['U'][fi::,k,:,:].squeeze()
                 v = nc.variables['V'][fi::,k,:,:].squeeze()
                 [nt, ny, nx] = u.shape
                 nc.close()
                 nc = Dataset (listadata02[n])
                 T = nc.variables['T'][fi::,k,:,:].squeeze()
                 #print ('nt = ', nt)
                 uc[ni:ni+nt,:,:] = u[:]
                 vc[ni:ni+nt,:,:] = v[:]
                 tc[ni:ni+nt,:,:] = T[:]
                 nc.close()
                 #print (num2date(ni,units=unitsd),ni) 
                 #print (num2date(ni+nt,units=unitsd),nt)
                 del u, v, T
                 ni = ni + nt
        except:
            print (syear[0])
    print ()
    print ('-----> high pass filter begins')

    [nt,ny,nx] = uc.shape
    print (n)
    #np.savetxt('prueba.dat', uc[:,0,0], fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
    # calculando los filtros para uc y vc del nivel seleccionado
    uf = np.zeros((nt,ny,nx))*np.nan
    vf = np.zeros((nt,ny,nx))*np.nan
    tf = np.zeros((nt,ny,nx))*np.nan
    print (uc[0,0,0],vc[0,0,0],tc[0,0,0])
    for j in range(0,ny):
        for i in range (0,nx):
            uf[:,j,i] = signal.filtfilt(B,A, uc[:,j,i])
            vf[:,j,i] = signal.filtfilt(B,A, vc[:,j,i])
            tf[:,j,i] = signal.filtfilt(B,A, tc[:,j,i])
    print (uf[0,0,0],vf[0,0,0],tf[0,0,0])
    uv = uf*vf
    vt = tf*vf
    ##%% esta parte es para la diferencial (inicia)
    print ()
    print ('-----> flux moment convergence begins for '+'level ' + str(nlev[k]) + ' mb')
    flux = uv * cos2
    # la derivada de la convergencia de flujo de momento
    dy = np.zeros(uv.shape,np.float)
    dy[:,0:-1,:] = ( np.diff(flux,axis=1) /np.diff(yd,axis=0) )
    # la derivada de la frontera
    dy[:,-1,:] = ( (flux[:,-1,:] - flux[:,-2,:]) /(yd[-1] - yd[-2]) )

    dy =  ( -1.0/(R*cos2) ) * dy

    print ()
    print ('-----> Done... derivative begins... ')
    print ()
    print ('-----> Now monthly means')
    a = datetime.date(year,1,1)
    idi = 0
    ti = []
    for t in range (0,nt,1):
        lm = monthrange(a.year,a.month)[1]
        if lm == 29:
            lm = 28
        idf = idi + lm
        #fdate = idf-10
        # guardando la derivada meridional
        # promedio zonal
        #var = np.nanmean(dy,axis=-1)
        # promedio mensual
        vardif = np.nanmean(dy[idi:idf,:,:], axis=0)
        # guardando la variable VU flujo de momento
        # promedio zonal
        #var = np.nanmean(uv,axis=-1)
        # promedio mensual
        varuv = np.nanmean(uv[idi:idf,:,:], axis=0)
        # guardando la variable VT flujo de calor
        # promedio mensual
        #var = np.nanmean(vt,axis=-1)
        # promedio mensual
        varvt = np.nanmean(vt[idi:idf,:,:], axis=0)
        FluxWind[nm,k,:,:] = vardif
        MomeWind[nm,k,:,:] = varuv
        Fluxtemp[nm,k,:,:] = varvt

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
print ()
#sio.savemat ('rediffFlux.mat', {'FluxWind':FluxWind,'lat':y,'lon':x,'lev':nlev})
os.system("rm diffFluxless.nc")
# abriendo el archivo netcdf
f = Dataset('diffFluxless.nc','w', format='NETCDF4') #'w' stands for write

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

EMFC = f.createVariable('EMFC', 'f4', ('time','lev','lat', 'lon')) #float
EMF = f.createVariable('EMF', 'f4', ('time','lev','lat', 'lon')) #float
MTF = f.createVariable('MTF', 'f4', ('time','lev','lat', 'lon')) #float

## long names
lon.long_name = "longitude"
lat.long_name = "latitude"
lev.long_name = "pressure_level"

EMFC.long_name = "Eddy Momentum Flux Convergence"
EMF.long_name = "Eddy Momentum Flux"
MTF.long_name = "Meridional temperature flux"

## units
lon.units  = "degrees_east"
lat.units  = "degrees_north"
lev.units  = "millibars"
time.units = unitsd
EMFC.units = "m/s^2"
EMFC.units = "m^2/s^2"
MTF.units  = "k*m/s"

## time calendar
time.calendar = "standard"

## asignando valores
lon[:]  = x[:]
lat[:]  = y[:]
lev[:]  = nlev[:]
time[:] = ti[:]
EMFC[:] = FluxWind[:]
EMF[:]  = MomeWind[:]
MTF[:]  = Fluxtemp[:]

f.history = "calculos con datos de experimento "+exp_case

f.close()

print ()
print ('Program done... you can check now the mat file with Eddy Momentum Flux Convergence [EMFC]')
toc = pytime.clock()
print
print ('time used in this process '+ str((toc -tic)/(60*60)) + ' hrs')
print ()
print ('Bye... ')
                         
