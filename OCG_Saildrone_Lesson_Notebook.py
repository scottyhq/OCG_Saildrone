# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # 2021 Oceanography Camp for Girls Saildrone Lesson
# Developed by Nancy Williams, Veronica Tamsitt, Nicola Guisewhite at University of South Florida College of Marine Science

# ## To Do List:
# * make a function for plotting instead of copy/pasting the map each time (or not, if we want to keep it simple)
# * explore https://jupytext.readthedocs.io/en/latest/ to iron out issues with version controlling Jupyter notebooks
# * add more markdown in the form of instructions, pictures, pulling variables out into their own cell so girls know where they can make changes to the code
# * in figure titles and filenames, change the variables from using the first four characters (currently var[:4]) to instead cutting off at the first underscore

# ## Data Sources:
# * Saildrone 1-minute physical and ADCP data available from: https://data.saildrone.com/data/sets/antarctica-circumnavigation-2019
# (login required, so cannot be accessed using an FTP. Will need to download ahead)
# * Saildrone hourly-ish CO2, pH data available from: https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0221912
# * Satellite Chlorophyll: https://neo.sci.gsfc.nasa.gov/view.php?datasetId=MY1DMW_CHLORA&year=2019
# * SSH: https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-sea-level-global?tab=overview 
# (login required for chla and SSH, download ahead of time. Can also be downloaded using motuclient, login also required https://github.com/clstoulouse/motu-client-python)

# Import the tools you need
import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature

# Set the paths
output_dir = 'Output/'
data_dir = 'Data/'

# Go and download the hourly Saildrone CO2 data and put it in the `Data/` folder
os.chdir(data_dir) # Change the directory to the `Data/` folder
os.getcwd() # Check that you're now in the `Data/` folder
# Curl downloads the data files directly from the web and shows you the status while it works. 
# `!` at the beginning of the line tells you that this command is a unix shell command (not python code)
# ! curl -o 32DB20190119_ASV_Saildrone1020_Antarctic_Jan2019_Aug2019.csv https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0221912/32DB20190119_ASV_Saildrone1020_Antarctic_Jan2019_Aug2019.csv
os.chdir("..") # Use ".." to move back up one directory now that we've imported the MLD climatology data

# Import the hourly Saildrone CO2 data file
Saildrone_CO2 = pd.read_csv(
    (data_dir + '32DB20190119_ASV_Saildrone1020_Antarctic_Jan2019_Aug2019.csv'),
    header=4,
    na_values=-999,
)
# Check that the Saildrone data was imported correctly
Saildrone_CO2

# Import the one-minute resolution Saildrone Physical data file
ds = xr.open_dataset(data_dir + 'saildrone-gen_5-antarctica_circumnavigation_2019-sd1020-20190119T040000-20190803T043000-1_minutes-v1.1620360815446.nc')
Saildrone_phys = ds.to_dataframe()
Saildrone_phys

# Import the Southern Ocean fronts for mapping
stf = pd.read_csv(data_dir + 'fronts/stf.txt', header=None, sep='\s+', 
                  na_values='%', names=['lon','lat'])
saf = pd.read_csv(data_dir + 'fronts/saf.txt', header=None, sep='\s+', 
                  na_values='%', names=['lon','lat'])
pf = pd.read_csv(data_dir + 'fronts/pf.txt', header=None, sep='\s+', 
                 na_values='%', names=['lon','lat'])
saccf = pd.read_csv(data_dir + 'fronts/saccf.txt', header=None, sep='\s+', 
                    na_values='%', names=['lon','lat'])
sbdy = pd.read_csv(data_dir + 'fronts/sbdy.txt', header=None, sep='\s+', 
                   na_values='%', names=['lon','lat'])

# +
# Plot the Saildrone track on a map

# Make the "bones" of the figure
plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.SouthPolarStereo())
ax.set_extent([-180, 180, -90, -30],ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND)
ax.add_feature(cartopy.feature.OCEAN, color='lightblue')
ax.gridlines()

# Compute a circle in axes coordinates, which we can use as a boundary
# for the map. We can pan/zoom as much as we like - the boundary will be
# permanently circular.
theta = np.linspace(0, 2 * np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

# Plot the ACC fronts in various colors
ax.set_boundary(circle, transform=ax.transAxes)
plt.plot(stf['lon'], stf['lat'], color='Red', transform=ccrs.PlateCarree(), 
         label = 'Subtropical Front')
plt.plot(saf['lon'], saf['lat'], color='Orange', transform=ccrs.PlateCarree(), 
         label = 'Subantarctic Front')
plt.plot(pf['lon'], pf['lat'], color='Yellow', transform=ccrs.PlateCarree(), 
         label = 'Polar Front')
plt.plot(saccf['lon'], saccf['lat'], color='Green', transform=ccrs.PlateCarree(), 
         label = 'Southern ACC Front')
plt.plot(sbdy['lon'], sbdy['lat'], color='Blue', transform=ccrs.PlateCarree(), 
         label = 'Southern Boundary of ACC')

# Plot the Saildrone in black dots
plt.scatter(Saildrone_phys.longitude, Saildrone_phys.latitude,
           transform=ccrs.PlateCarree(), c='black', s=3, label='Saildrone', zorder=1000)

# Turn on the legend
plt.legend()

# Save the figure in the output folder
plt.title('2019 Saildrone Antarctic Circumnavigation Track')
plt.savefig(output_dir + 'SaildroneMap' + '.jpg') # Changing the suffix will change the format
plt.show()

# +
# Now plot some variable "var" on the map with colored dots
var = 'TEMP_CTD_RBR_MEAN'

# Make the "bones" of the figure
plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.SouthPolarStereo())
ax.set_extent([-180, 180, -90, -30],ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND)
ax.add_feature(cartopy.feature.OCEAN, color='lightblue')
ax.gridlines()

# Compute a circle in axes coordinates, which we can use as a boundary
# for the map. We can pan/zoom as much as we like - the boundary will be
# permanently circular.
theta = np.linspace(0, 2 * np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

# Plot the ACC fronts in various colors
ax.set_boundary(circle, transform=ax.transAxes)
plt.plot(stf['lon'], stf['lat'], color='Red', transform=ccrs.PlateCarree(), 
         label = 'Subtropical Front')
plt.plot(saf['lon'], saf['lat'], color='Orange', transform=ccrs.PlateCarree(), 
         label = 'Subantarctic Front')
plt.plot(pf['lon'], pf['lat'], color='Yellow', transform=ccrs.PlateCarree(), 
         label = 'Polar Front')
plt.plot(saccf['lon'], saccf['lat'], color='Green', transform=ccrs.PlateCarree(), 
         label = 'Southern ACC Front')
plt.plot(sbdy['lon'], sbdy['lat'], color='Blue', transform=ccrs.PlateCarree(), 
         label = 'Southern Boundary of ACC')

# Plot the Saildrone in black dots
plt.scatter(Saildrone_phys.longitude, Saildrone_phys.latitude, 
            c=Saildrone_phys[var], cmap='bwr',
            transform=ccrs.PlateCarree(), s=5, zorder=1000)

# Turn on the legend
plt.legend()
cb1 = plt.colorbar()

# Save the figure in the output folder
plt.title('2019 Saildrone ' + var[:4])
plt.savefig(output_dir + var[:4] + 'SaildroneMap' + '.jpg') # Changing the suffix will change the format
plt.show()
# -

# Now let's make a map of some satellite data to show the Saildrone crossing an ocean eddy.
#
# First we need to load in a single daily satellite sea surface height data file from Feb 10th 2019, the day the Saildrone crossed a large eddy.

satellite_ssh = xr.open_dataset(data_dir + 'ssh_2019_02_10.nc')

# Now plot the Saildrone path on a map of sea surface height for a region surrounding the Saildrone on Feb 10th

#set plot parameters (contour levels, colormap etc)
levels_1 = np.arange(-1.2,0.8,0.1) #contour levels
cmap_1 = 'viridis' #contour map colormap
c1 = 'black' #Saildrone track color

# +
#finding position of Saildrone on Feb 10
time_index = np.argwhere(Saildrone_phys.time.values==np.datetime64('2019-02-10'))[0]
tlon = Saildrone_phys.longitude.values[time_index]
tlat = Saildrone_phys.latitude.values[time_index]

#make a contour plot of satellite ssh
xr.plot.contourf(satellite_ssh.adt[0,:,:],levels=levels_1,cmap=cmap_1,size=8,aspect=2)
xr.plot.contour(satellite_ssh.adt[0,:,:],levels=levels_1,colors='k',linewidths=0.75)
plt.xlim(tlon+360-5,tlon+360+5)
plt.ylim(tlat-5,tlat+5)

#add Saildrone track
plt.scatter(Saildrone_phys.longitude+360, Saildrone_phys.latitude, c=c1, s=3, label='Saildrone', zorder=1000)
plt.legend()

#give the plot a title and save figure in the output folder
plt.title('Saildrone path across an eddy on Feb 10th')
plt.savefig(output_dir + 'Sea_surface_height_Saildrone_Feb10' + '.jpg')

# -

# Ocean eddies can be identified by closed rings of constant absolute dynamic topography (this is the anomaly in sea surface height from average sea level in meters, which represents changes in pressure). You can see the Saildrone's path crossing near the center of an eddy.
#
# We can do the same thing with satellite chlorophyll-a data. The chlorophyll-a data gives an approximate estimate of the relative phytoplankton biomass (in units of mg/m<sup>3</sup>) at the sea surface in different locations. 

#set plot parameters (contour levels, colormap etc)
levels_1 = np.arange(0,1.0,0.01) #contour levels
cmap_1 = 'YlGnBu' #contour map colormap
c1 = 'black' #Saildrone track color

# +
#load satellite chl-a data file
satellite_chla = xr.open_dataset(data_dir + 'A2019041.L3m_DAY_CHL_chlor_a_4km.nc')

#make a contour plot of chl-a data 
xr.plot.contourf(satellite_chla.chlor_a,levels=levels_1,cmap=cmap_1,size=8,aspect=2)
plt.xlim(tlon-30,tlon+30)
plt.ylim(tlat-20,tlat+20)

#add Saildrone track
plt.scatter(Saildrone_phys.longitude, Saildrone_phys.latitude, c=c1, s=3, label='Saildrone', zorder=1000)
plt.legend()

#save figure
plt.title('Saildrone path and chlorophyll-a concentration')
plt.savefig(output_dir + 'Sea_surface_chlorophylla_Saildrone_Feb10' + '.jpg')
# -

# Now we can add the Saildrone data observations on the map to start to see if there is a relationship between the satellite observations and what the Saildrone measured directly. Note that the Saildrone took a few days to cross this region, while the satellite data shown here is a snapshot for a single day, so it can be tricky to compare the two types of data because the Saildrone is moving in space AND time.

#choose which variable to plot
var = 'TEMP_CTD_RBR_MEAN'
#set minimum and maximum colorbar limits
v_min = 6
v_max = 12
#choose colormap for map
cmap_1 = 'viridis'
#choose colormap for Saildrone variable
cmap_2 = 'RdBu_r'

# +
#make a contour plot of satellite ssh
xr.plot.contourf(satellite_ssh.adt[0,:,:],levels=np.arange(-1.2,0.8,0.1),cmap=cmap_1,size=8,aspect=2)
xr.plot.contour(satellite_ssh.adt[0,:,:],levels=np.arange(-1.2,0.8,0.1),colors='black',linewidths=0.75)
plt.xlim(tlon+360-5,tlon+360+5)
plt.ylim(tlat-5,tlat+5)

#add Saildrone data scattered on top
plt.scatter(Saildrone_phys.longitude+360, Saildrone_phys.latitude, c=Saildrone_phys[var], s=15, cmap = cmap_2,
            vmin=v_min,vmax=v_max,label='Saildrone', zorder=1000)
plt.legend()
plt.colorbar()

#add title and save figure
plt.title('Saildrone ' + var[:4] + ' across an eddy on Feb 10th')
plt.savefig(output_dir + 'Sea_surface_height_Saildrone_' + var[:4] + '_Feb10' + '.jpg')
# -

# Now that we've plotted some Saildrone and satellite data on maps to see how different ocean variables are related, there are other ways we can look at the relationship between variables. 
#
# This includes scatter plots, which is a useful way to compare data from two variables collected at the same time and location to look for a relationship. In our case, we can compare two different variables collected simulatneously by the Saildrone. 

# +
#choose two variables from the Saildrone to compare
var1 = 'TEMP_CTD_RBR_MEAN'
var2 = 'O2_CONC_RBR_MEAN'

#choose lower and upper limits for the two variables for plotting
var1_min = -2
var1_max = 18
var2_min = 230
var2_max = 340

#choose color for scatterplot
c1 = 'b'

# +
#create scatter plot
plt.figure(figsize=(12,8))
plt.scatter(Saildrone_phys[var1], Saildrone_phys[var2], c=c1, s=10)
plt.xlim(var1_min,var1_max)
plt.ylim(var2_min,var2_max)
plt.xlabel(var1)
plt.ylabel(var2)
plt.grid()

#add title and save figure
plt.title('Saildrone '+ var1[:4] + ' vs ' + var2[:4])
plt.savefig(output_dir + 'Saildrone_' + var1[:4] + '_vs_' + var2[:4] + '.jpg')
# -

# If we want to look at the relationship between more than two variables, one way we can look at this is by using a third variable to change the color of the scatter plot points.
#
# In this example, the scatter plot shows the same variable 1 and 2 from the Saildrone on the x and y axes as above, but we can choose a third Saildrone variable as the color of the scatter plot points.

# +
#choose a third variable 
var3 = 'latitude'

#set lower and upper limits of variable 3 for plotting
var3_min = -65
var3_max = -40

# +
#create scatter plot
plt.figure(figsize=(12,8))
plt.scatter(Saildrone_phys[var1], Saildrone_phys[var2], c=Saildrone_phys[var3], 
            s=10, vmin = var3_min, vmax = var3_max)
plt.xlim(var1_min,var1_max)
plt.ylim(var2_min,var2_max)
plt.xlabel(var1)
plt.ylabel(var2)
plt.grid()
cbar = plt.colorbar()
cbar.set_label(var3)

#add title and save figure
plt.title('Saildrone '+ var1[:4] + ' vs ' + var2[:4])
plt.savefig(output_dir + 'Saildrone_' + var1[:4] + '_vs_' + var2[:4] + '_vs_' + var3[:4] + '.jpg')
# -

# Plot time series of wind and pressure data

#set plot parameters
var1 = 'UWND_MEAN'
var2 = 'BARO_PRES_MEAN'

# +
#plot time series
plt.figure(figsize=(12,5))
ax1 = plt.subplot(211)
ax1.plot(Saildrone_phys.time,Saildrone_phys[var1])
plt.xlim(Saildrone_phys.time.values[0],Saildrone_phys.time.values[-1])

ax2 = plt.subplot(212)
ax2.plot(Saildrone_phys.time,Saildrone_phys[var2])
plt.xlim(Saildrone_phys.time.values[0],Saildrone_phys.time.values[-1])
# -

