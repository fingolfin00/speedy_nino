#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 10:41:30 2022

@author: Paolo
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy

proj=ccrs.PlateCarree(central_longitude=180)
proj0=ccrs.PlateCarree(0)

seas='DJF'
DIR='./Data/' #path to directory with NetCDF files 
FILEclim_djf=DIR+seas+'mean_clim.nc' #Full path with file name
FILEnino_djf=DIR+seas+'mean_nino.nc' #Full path with file name

seas='JJA'
DIR='./Data/' #path to directory with NetCDF files 
FILEclim_jja=DIR+seas+'mean_clim.nc' #Full path with file name
FILEnino_jja=DIR+seas+'mean_nino.nc' #Full path with file name

ncin=Dataset(FILEclim_djf,'r')
chi_djf=ncin.variables['chi'][:]
psi_djf=ncin.variables['psi'][:]
u_djf=ncin.variables['u'][:]
v_djf=ncin.variables['v'][:]
tp_djf=ncin.variables['tp'][:]
temp0_djf=ncin.variables['temp0'][:]
omega_djf=ncin.variables['omega'][:]
lon=ncin.variables['lon'][:]
lat=ncin.variables['lat'][:]
lev=ncin.variables['lev'][:]
ncin.close()

ncin=Dataset(FILEclim_jja,'r')
chi_jja=ncin.variables['chi'][:]
psi_jja=ncin.variables['psi'][:]
u_jja=ncin.variables['u'][:]
v_jja=ncin.variables['v'][:]
tp_jja=ncin.variables['tp'][:]
temp0_jja=ncin.variables['temp0'][:]
omega_jja=ncin.variables['omega'][:]
ncin.close()

ncin_nino=Dataset(FILEnino_djf,'r')
chi_nino_djf=ncin_nino.variables['chi'][:]
psi_nino_djf=ncin_nino.variables['psi'][:]
u_nino_djf=ncin_nino.variables['u'][:]
v_nino_djf=ncin.variables['v'][:]
tp_nino_djf=ncin_nino.variables['tp'][:]
temp0_nino_djf=ncin_nino.variables['temp0'][:]
omega_nino_djf=ncin_nino.variables['omega'][:]
ncin_nino.close()

ncin_nino=Dataset(FILEnino_jja,'r')
chi_nino_jja=ncin_nino.variables['chi'][:]
psi_nino_jja=ncin_nino.variables['psi'][:]
u_nino_jja=ncin_nino.variables['u'][:]
v_nino_jja=ncin.variables['v'][:]
tp_nino_jja=ncin_nino.variables['tp'][:]
temp0_nino_jja=ncin_nino.variables['temp0'][:]
omega_nino_jja=ncin_nino.variables['omega'][:]
ncin_nino.close()

#South America
las,lan=-25,-5
low,loe=290,320
#Europe
low,loe=0,40
las,lan=55,75
#Central America
low,loe=240,270
las,lan=20,40

def generate_zonal_fields(field_jja, field_djf, field_nino_jja, field_nino_djf):
    field_zonal = np.nanmean((field_jja+field_djf)/2, axis=3)
    field_zonal_mean = np.mean(field_zonal, axis=0)

    field_zonal_nino = np.nanmean((field_nino_jja+field_nino_djf)/2, axis=3)
    field_zonal_nino_mean = np.mean(field_zonal_nino, axis=0)

    return field_zonal_mean, field_zonal_nino_mean

u_zonal_mean, u_zonal_nino_mean = generate_zonal_fields(u_jja, u_djf, u_nino_jja, u_nino_djf)
v_zonal_mean, v_zonal_nino_mean = generate_zonal_fields(v_jja, v_djf, v_nino_jja, v_nino_djf)
omega_zonal_mean, omega_zonal_nino_mean = generate_zonal_fields(omega_jja, omega_djf, omega_nino_jja, omega_nino_djf)

x_lat, y_lev = np.meshgrid(lat, lev)

fig,ax=plt.subplots(1,1)

c = ax.contour(x_lat,y_lev,u_zonal_mean,colors='k',linewidths=1)
cf = ax.contourf(x_lat, y_lev, u_zonal_nino_mean - u_zonal_mean,
              cmap='seismic', levels=np.linspace(-1.3,1.3,12))
ax.clabel(c, inline=True, fontsize=10)

#c=ax[1].contourf(x_lat,y_lev,u_zonal_nino_mean,cmap='seismic')
plt.colorbar(cf, ax=ax, orientation='horizontal')
plt.gca().invert_yaxis()

plt.savefig("u_zonal_anomaly.png")

fig,ax=plt.subplots(1,1)

c = ax.contour(x_lat,y_lev,v_zonal_mean,colors='k',linewidths=1)
cf = ax.contourf(x_lat, y_lev, v_zonal_nino_mean - v_zonal_mean,
              cmap='seismic', levels=np.linspace(-0.15,0.15,12))
ax.clabel(c, inline=True, fontsize=10)

#c=ax[1].contourf(x_lat,y_lev,u_zonal_nino_mean,cmap='seismic')
plt.colorbar(cf, ax=ax, orientation='horizontal')
plt.gca().invert_yaxis()

plt.savefig("v_zonal_anomaly.png")

fig,ax=plt.subplots(1,1)

c = ax.contour(x_lat,y_lev,omega_zonal_mean,colors='k',linewidths=1, levels=np.linspace(-0.02,0.02,8))
cf = ax.contourf(x_lat, y_lev, omega_zonal_nino_mean - omega_zonal_mean,
              cmap='seismic', levels=np.linspace(-0.0025,0.0025,12))
ax.clabel(c, inline=True, fontsize=10)

#c=ax[1].contourf(x_lat,y_lev,u_zonal_nino_mean,cmap='seismic')
plt.colorbar(cf, ax=ax, orientation='horizontal')
plt.gca().invert_yaxis()

plt.savefig("omega_zonal_anomaly.png")