# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 13:31:20 2018

@author: Jose Luis Rodriguez-Solis
         Oceanografia Fisica
         CICESE

LAT POSITION
Calcula los valores de posicion latitudinal  basados en :

Archer, C. L., & Caldeira, K. (2008). Historical trends in the jet streams. 
Geophysical Research Letters, 35(8), 1–6. doi.org/10.1029/2008GL033614


"""

from netCDF4 import Dataset, date2num, num2date
import os
import scipy.io as sio
import numpy as np
import atmos
import datetime
from calendar import monthrange, isleap
from dateutil.relativedelta import relativedelta
import glob

#run = 'lamip'

exp_case='mamip'

f01=glob.glob('/scratch/201803026n-2/cesm/var/'+exp_case+'/var00/*.nc')
f01.sort()
f02=glob.glob('/scratch/201803026n-2/cesm/var/'+exp_case+'/var01/*.nc')
f02.sort()
f03=glob.glob('/scratch/201803026n-2/cesm/var/'+exp_case+'/var02/*.nc')
f03.sort()

# nivel superior
pi = 400 
# nivel infererior
ps = 100
# definiendo las constantes
Rd = 287.04 # kg * m^2 * s^-2
g  = 9.81   # m * s^-2

# donde estan los archivos de Era Interim
#path = '/HOME/users/mgross/cesm/FAMIP2x2_more/run/'
#f = os.listdir(path)
#f.sort()


nc = Dataset (f01[0])
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
lev = nc.variables['lev'][:]
timed     = nc.variables['time'][:]
unitsd    = nc.variables['time'].units
[nt, nz, ny, nx] = nc.variables['U'].shape
nc.close()
# buscando indices de altura de 400 y 100
zl  = np.where((lev>= 100.0) & (lev<=  400.5)); zl   = zl[0]
zs = zl[0]
zi = zl[-1]

# length period
dt = (datetime.date(2005, 12, 31) - datetime.date(1981, 1, 1)).days
# leap years count
yleap = 0
for i in range(1981,2006):
    if isleap(i):
        yleap = yleap + 1
dt = dt - yleap + 1

year = datetime.datetime(1981,1,1).year
si   = datetime.datetime(1980,12,1)
# Malla para graficar y encontrar latitudes
x,y=np.meshgrid(lon,lat)
# buscando latitudes entre 70N y 15N
yl  = np.where((lat>= 15.0) & (lat<=  70.5)); yl   = yl[0]
ys = yl[0]
yi = yl[-1]

xl  = np.where((lon>= 229.5) & (lon<=  250.5)); xl   = xl[0]
xs = xl[0]
xi = xl[-1]

#print (ys,yi,yl)
Dz = []
for i in range(0, len(lev)-1):    
    Dz.append((lev[i] + lev[i+1]) / 2.0)
Dz = np.asarray(Dz)

Dzd = []
for i in range(0, len(Dz)-1):    
    Dzd.append(Dz[i+1]/Dz[i])
Dzd = np.asarray(Dzd)

del Dz

zl = np.zeros((len(zl),len(lat),len(lon)))*np.nan
for i in range(0, len(lev[zs:zi+1])):    
    zl[i,:,:] = lev[zs+i]

Dz = (Rd/g)*np.log(Dzd[zs-1:zi])

dz = np.zeros((len(Dz),len(lat),len(lon)))*np.nan
for i in range(0, len(Dz)):
    dz[i,:,:] = Dz[i]

#sio.savemat ('dz.mat', {'dz':dz,'zl':zl})
print (dz.shape)

SLat   = np.zeros((12*25))

ws = np.zeros((dt,ny,nx))* np.nan
Ps = np.zeros((dt,ny,nx))* np.nan
ls = np.zeros((dt))* np.nan
#ls = []
print (ws.shape)
num = len(f01)
ni = 0
print (xs,xi)
for n in range(0,len(f01)):
    nc = Dataset (f01[n])
    timed = nc.variables['time'][:]
    syear = num2date(timed,units=unitsd)
   
    fi = np.where(syear>si)
    try:
        if syear[fi[0][0]] > si:
            fi = fi[0][0]
            print (syear[fi])  
            nc = Dataset (f01[n])    
            u   = nc.variables['U'][fi::,zs:zi+1,:,:].squeeze()
            v   = nc.variables['V'][fi::,zs:zi+1,:,:].squeeze()
            us = u.copy()**2
            vs = v.copy()**2
            #print (u[0,0,0,0],v[0,0,0,0])
            del u, v
            nc.close()
            nc = Dataset (f02[n])
            T   = nc.variables['T'][fi::,zs:zi+1,:,:].squeeze()
            nc.close()
            nc = Dataset (f03[n])    
            qv  = nc.variables['Q'][fi::,zs:zi+1,:,:].squeeze()
            nc.close()    
            #print (qv[0,0,0,0],T[0,0,0,0])
            [nt, nz, ny, nx] = us.shape
            #print ('ni, nt = ',ni, nt)
            # buscando indices de niveles superior e inferior
            ms = np.sqrt(us + vs)
            ts =  T.copy()
            del T
            # calculando rho
            Tv = (1 + 0.608*qv)*ts
            rho = atmos.calculate('rho',p = zl*100.0, Tv = Tv)
            # calculando el espesor de la capa en metros
            dl = dz * Tv
            rho = rho * dl
            #calculando latitud
            del qv
            ls_top = np.nansum ( rho * ms, axis=1) * y
            #print (ls_top.shape)
            ls_top = np.nansum ( ls_top[:,ys:yi+1,xs:xi+1], axis = 1)
            #print (ls_top.shape)
            ls_top = np.nanmean ( ls_top, axis=1)
            #print (ls_top.shape)

            ls_inf = np.nansum ( rho * ms, axis=1)
            ls_inf = np.nansum ( ls_inf[:,ys:yi+1,xs:xi+1], axis = 1)
            ls_inf = np.nanmean ( ls_inf, axis=1)
            print (ls_inf.shape,ls_top.shape) 
            
            #lstot = np.nanmean(ls_top/ls_inf) 
            #print (lstot)
            ls[ni:ni+nt] = ( ls_top/ls_inf )
            #ls.append(lstot)
            ni = ni + nt
    except:
        print ('no entra') 
          
print ('finaliza ciclo de calculos')
print ('')
print ('comienza ciclo para medias')
#%% inicia promedios 

ls = np.asarray(ls)

#%% PROMEDIOS MENSUALES
a = datetime.date(year,1,1)
idi = 0
ti = []
nm = 0
for t in range (0,len(ws)+1,1):
    lm = monthrange(a.year,a.month)[1]
    if lm == 29:
       lm = 28
    idf = idi + lm
    # promedio mensual
    varvt = np.nanmean(ls[idi:idf], axis=0)
    SLat[nm]   = varvt * 1.0
    print ('>> month number :' + str(nm) + ' ' + a.strftime('%Y-%m'))
    nm = nm+1
    b  = datetime.datetime(a.year,a.month,a.day)
    b = date2num(b,unitsd,calendar='standard')
    ti.append(b)
    a = a + relativedelta(months=1)
    print (t,idi,idf,nt)
    idi = idf*1
    if a.year == 2006 and a.month == 1:
        break

#%% serie de tendencia de Presion 
ti  = np.asarray(ti)
tip = np.arange(1981,2005+1,1)

sio.savemat ('Lat_'+exp_case+'09x125.mat', {'SLat':SLat,'t01':ti,'t02':tip,'dz':dz})

print ('Fin del programa')
