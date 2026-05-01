
## for plotting SSH / ocean color underway with ship track 

import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from datetime import datetime
from os import listdir
import glob
import os
import pandas as pd
import imageio
import scipy
from netCDF4 import Dataset
import datetime as dt
import cartopy.crs as ccrs

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy
import earthaccess
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import xarray as xr
import os
import netCDF4 
import fnmatch

from matplotlib.animation import FuncAnimation, PillowWriter  
import matplotlib.colors as mcolors

path = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/Kahru_composites/'
path = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/Kahru_composites/'
products = ['Chl5day1km', 'Chl15day1km', 'ChlMo500m', 'ChlMo1km', 'Chlday1km']
files = os.listdir(path+products[0]+'/')

import pyhdf.SD as phs
#from pyhdf.SD import SD, SDC

# product of interest: 15 day composite from Feb 15 (DOY = 46) to March 1 (DOY = 60) 

path = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/Kahru_composites/'
fn = 'C20250462025060_chl_comp.hdf'

f = phs.SD(path+fn, phs.SDC.READ)

print(f.datasets()) 
chl = f.select('Chl') 
data = chl[:,:].astype(float)
data.shape


# also check out the mapping file 
fn = 'cal_aco_3840_Latitude_Longitude.hdf'
hmap = path+fn

f = phs.SD(hmap, phs.SDC.READ)

print(f.datasets()) 
variables = f.datasets()
print(variables)
lat_name = list(variables.keys())[0]
lon_name = list(variables.keys())[1]

lons = f.select(lon_name) 
lats = f.select(lat_name) 
lon = lons[:,:]; lat = lats[:,:]

# process according to mati's explanation: https://spg-satdata.ucsd.edu/Readme.htm
# when reading with Matlab the unsigned byte variable is reported as signed byte (int8) and values over 128 become negative. A simple fix is to add 256 if the signed pixel value is negative.

data[data<0] += 256


plt.figure()
plt.pcolormesh(lon,lat,data)
plt.colorbar()
plt.show()


# Pixel values 0 and 255 (and the corresponding scaled values) are considered invalid and must be excluded in any statistics.

data[data==0] = np.nan
data[data==255] = np.nan

plt.figure()
plt.pcolormesh(lon,lat,data)
plt.colorbar()
plt.show()

# sst & chl in hdf use 1 byte per pixel with specific scaling. Linear scaling is used for SST and logarithmic scaling for Chl. The scaling equations using pixel value (PV) as unsigned byte (from 0 to 255) are:
# SST (deg C) = 0.15 * PV - 3.0
# Chl (mg m-3) = 10^(0.015 * PV - 2.0), i.e. 10 to the power of 0.015 * PV - 2.0.

power = 0.015*data - 2.0
chl = np.power(10,power)
logchl = np.log(chl)


plt.figure()
plt.pcolormesh(lon,lat,logchl)
plt.colorbar()
plt.show()

# actually - use each station instead 

path = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/SR2503_scienceshare/CTD/Survey_Stations/'
files = os.listdir(path)


def parse_cnv_file(filename):
    """Parse a Sea-Bird .cnv file into metadata and data."""
    metadata = {}
    data = []

    with open(filename, 'r') as file:
        for line in file:
            # Check for metadata (lines starting with '*')
            if line.startswith('*'):
                if '=' in line:
                    # Metadata key-value pair
                    key, value = line[1:].strip().split('=', 1)
                    metadata[key.strip()] = value.strip()
                else:
                    # Metadata comments without key-value
                    key = line[1:].strip()
                    metadata[key] = None
            elif line.startswith('#'):
                # Column headers start with '#'
                if 'name' in line:
                    _, col_index, col_name = line.strip().split(' ', 2)
                    metadata[f"Column_{col_index.strip()}"] = col_name.strip()
            else:
                # Data starts after metadata and column definitions
                data.append(line.strip())

    # Convert data into a DataFrame
    if data:
        data = [list(map(float, row.split())) for row in data if row]
        df = pd.DataFrame(data, columns=[metadata.get(f"Column_{i}", f"Column_{i}") for i in range(len(data[0]))])
    else:
        df = pd.DataFrame()

    return metadata, df

def convert_to_decimal_degrees(coord_str):
    """
    Converts strings like '33 49.38 N' or '118 37.74 W' to decimal degrees.
    """
    parts = coord_str.split()
    degrees = float(parts[0])
    minutes = float(parts[1])
    direction = parts[2]

    decimal = degrees + minutes / 60.0
    if direction in ['S', 'W']:
        decimal = -decimal
    return decimal
    
def extract_lat_lon_time(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Ensure file has at least 13 lines
    if len(lines) < 13:
        raise ValueError("File does not contain enough lines.")

    # Strip '*' and whitespace, then split at '=' to isolate values
    lat_line = lines[9].strip().lstrip('*').strip()
    lon_line = lines[10].strip().lstrip('*').strip()
    time_line = lines[11].strip().lstrip('*').strip()

    lat_str = lat_line.split('=')[1].strip()
    lon_str = lon_line.split('=')[1].strip()
    utc_str = time_line.split('=')[1].strip()

    # Convert to desired formats
    latitude = convert_to_decimal_degrees(lat_str)
    longitude = convert_to_decimal_degrees(lon_str)
    utc_time = datetime.strptime(utc_str, "%b %d %Y %H:%M:%S")

    return latitude, longitude, utc_time

# find cast files 

matched_files = [f for f in files if fnmatch.fnmatch(f, '*_down.cnv')]

# loop through casts 1-13

stations = np.arange(1,14) 
num_stations = len(stations)

slons = []
slats = []
sdates = []


num_variables = 17
variable_order = ['Temp','Cond','Pressure','Temp2','Cond2','Transmissometer','Fluorometer','Alimeter','PAR','Oxygen','UPV','SPAR','Salinity','Salinity2','Lon','Lat','Flag']
# longest (?) shape of p/t/c object: 212103
deepest = 222873
#stations = [1,5,9,11,7,2,6,10,12,8,4,3];
#cardinallocs = np.array([0,9,5,1,3,7,10,6,2,4,8,12,11]) 
#cardinal = cardinallocs[1:]

survey_data = np.zeros((num_stations,num_variables,deepest))*np.nan


for ii,stat in enumerate(stations): 
    # save lat/lon
    filename = f'SR2503_cast_{stat:02d}_updown.cnv'
    [this_lat, this_lon, this_time] = extract_lat_lon_time(path+filename)
    slons.append(this_lon)
    slats.append(this_lat)
    sdates.append(this_time)

    # save ctd/bio variables
    metadata, data = parse_cnv_file(path+filename)

    # save data with index as each station #
    zlen = len(data['Column_2'].values)
    
    temp1 = data['Column_0'].values
    survey_data[ii,0,0:zlen] = temp1
    
    cond1 = data['Column_1'].values
    survey_data[ii,1,0:zlen] = cond1

    p = data['Column_2'].values
    survey_data[ii,2,0:zlen] = p

    temp2 = data['Column_3'].values
    survey_data[ii,3,0:zlen] = temp2

    cond2 = data['Column_4'].values
    survey_data[ii,4,0:zlen] = cond2
    
    trans = data['Column_5'].values
    survey_data[ii,5,0:zlen] = trans
    
    fluoro = data['Column_6'].values
    survey_data[ii,6,0:zlen] = fluoro
    
    altimeter = data['Column_7'].values
    survey_data[ii,7,0:zlen] = altimeter

    PAR = data['Column_8'].values
    survey_data[ii,8,0:zlen] = PAR
    
    oxygen = data['Column_9'].values
    survey_data[ii,9,0:zlen] = oxygen
    
    upv = data['Column_10'].values
    survey_data[ii,10,0:zlen] = upv
    
    SPAR = data['Column_11'].values
    survey_data[ii,11,0:zlen] = SPAR
    
    sal1 = data['Column_12'].values
    survey_data[ii,12,0:zlen] = sal1
    
    sal2 = data['Column_13'].values
    survey_data[ii,13,0:zlen] = sal2

    print(filename)
  


NN = len(sdates)



station_names = [1,5,9,11,7,2,6,10,12,8,4,3];

# now limit to region of interest: 
import matplotlib as mpl

title = 'Santa Lucia 15-day composite CHL \n Feb 15-March 1'

# set bounds
xmin, xmax = -124, -117
ymin, ymax = 32.5, 36
sdlonlat = [-117.24556, 32.928333]
chl_extent = [-4,3]

fig, ax = plt.subplots(figsize=(18, 8))
mesh = ax.pcolormesh(lon, lat, logchl, vmin=chl_extent[0], vmax=chl_extent[1]) #, cmap="jet"
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_title(title)

# log scale plot: 
# change mins: 
#chl_extent = [1e-2, 3]
#mesh = ax.pcolormesh(lon,lat,chl_data[0],cmap='jet',transform=crs_data,norm=mcolors.LogNorm(vmin=chl_extent[0],vmax=chl_extent[1]) )
cbar = plt.colorbar(mesh, ax=ax, orientation='vertical', pad=0.05)
cbar.set_label(r'Log Chlorophyll Concentration [mg/$m^3$]')


# add ship trajectory
numstops = np.arange(0,NN)
j = ax.plot(slons,slats,linewidth=1,label='Ship Track',color='k')
traj = ax.scatter(slons,slats,c=numstops,s=40,linewidth=1,cmap='jet')

j1 = ax.scatter(slons[0],slats[0],c='tab:red',marker='*',s=60)
ax.annotate('Sampling\nStart',(slons[0]-0.075,slats[0]-0.25),color='tab:red',fontsize=8,fontweight='bold')


# trajectory colorbar
cax1 = ax.inset_axes([0.2, 1.04, 0.8, 0.04]) # [left, bottom, width, height]
trajcbar = fig.colorbar(traj, cax=cax1,orientation='horizontal')
cax1.tick_params(axis='both',which='both',right=True,left=False,labelright=True)
cax1.xaxis.set_ticks_position('top'); cax1.xaxis.set_label_position('top') 
dates = [dt.datetime.strftime(d,'%m-%d') for d in sdates]
cax1.xaxis.set_ticks(numstops,labels=dates)
trajcbar.set_label('Dates along track',labelpad=5,fontsize=10)

# reference SD 
ax.scatter(sdlonlat[0], sdlonlat[1], color="red", marker="*", label="SD")

plt.show()

'''
fig,ax = plt.subplots(figsize=(12,6))

ymin = 4200 
base_cmap = plt.get_cmap('jet')  
colors = base_cmap(np.linspace(0, 1, NN))
a = 1
ymin = 200

for ii,istat in enumerate(stations[1:]): 
    p_atstat = survey_data[ii,2,:]
    c_atstat = survey_data[ii,6,:]
    ax.plot(c_atstat,p_atstat,label=f'Station {istat}',alpha=a,color=colors[ii])
    a = a - 0.05
    ax.set_ylim([ymin,0]); ax.set_xlabel('Fluorescence [spt]')
    ax.legend()
    ax.grid(True); plt.tick_params(axis='y',which='both',left=False,right=False, labelleft=False)

'''

# chl / fluoro 
fig,axes = plt.subplots(1,12,figsize=(20,4))

ymin = 100
var_cmap = mpl.cm.jet
var_colors = var_cmap(np.linspace(0,1,13))
a = 1 
for ii,istat in enumerate(stations[1:]): 
    ax = axes[ii]
    p_atstat = survey_data[ii+1,2,:]
    c_atstat = survey_data[ii+1,6,:]
    ax.plot(c_atstat,p_atstat,label=f'Station {istat}',color=var_colors[ii+1])
    a = a - 0.025
    ax.set_ylim([ymin,0]);
    ax.set_xlim([0,15])
    if ii == 5: 
        ax.set_xlabel('Fluorescence [spt]')
    if ii == 0 : 
        ax.tick_params(axis='y',which='both',left=True,right=False, labelleft=True)
        ax.set_ylabel('Depth [m]')
    else: 
        ax.tick_params(axis='y',which='both',left=False,right=False, labelleft=False)
    ax.grid(True); 
    ax.set_title(f'Station #{station_names[ii]}')
    
    # change color of all 4 borders (spines)
    for spine in ax.spines.values():
        spine.set_edgecolor(var_colors[ii+1])  # or any other color, e.g., '#00FF00'

plt.show()

