#!/usr/bin/env python
# coding: utf-8

# # In Class Lab 3
# 

# In[1]:


#In Class Lab 3 Template
# G Besla ASTR 400B

# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')


# The Figure illustrates the color magnitude diagram (CMD) for the Carina Dwarf along with the interpreted 
# star formation history from isochrone fitting to the CMD.
# The image is from Tolstoy+2009 ARA&A 47 review paper about dwarf galaxies
# 
# ![Iso](./Lab3_Isochrones.png)
# 

# # This Lab:
# 
# Modify the template file of your choice to plot isochrones that correspond to the inferred star formation episodes (right panel of Figure 1) to recreate the dominant features of the CMD of Carina (left panel of Figure 1). 

# In[2]:


# Some Notes about the Isochrone Data
# DATA From   
# http://stellar.dartmouth.edu/models/isolf_new.html
# files have been modified from download.  
# ( M/Mo --> M;   Log L/Lo --> L)
# removed #'s from all lines except column heading
# NOTE SETTINGS USED:  
# Y = 0.245 default   [Fe/H] = -2.0  alpha/Fe = -0.2
# These could all be changed and it would generate 
# a different isochrone


# In[3]:


# Filename for data with Isochrone fit for 1 Gyr
# These files are located in the folder IsochroneData
filename1="./IsochroneData/Isochrone1.txt"


# In[4]:


# READ IN DATA
# "dtype=None" means line is split using white spaces
# "skip_header=8"  skipping the first 8 lines 
# the flag "names=True" creates arrays to store the date
#       with the column headers given in line 8 

# Read in data for an isochrone corresponding to 1 Gyr
data1 = np.genfromtxt(filename1,dtype=None,
                      names=True,skip_header=8)


# In[11]:


#major peak
filename11="./IsochroneData/Isochrone11.txt"
filename10="./IsochroneData/Isochrone10.txt"

data11 = np.genfromtxt(filename11,dtype=None,
                      names=True,skip_header=8)
data10 = np.genfromtxt(filename10,dtype=None,
                      names=True,skip_header=8)


# In[20]:


#next peak
filename6="./IsochroneData/Isochrone6.txt"
filename7="./IsochroneData/Isochrone7.txt"

data6 = np.genfromtxt(filename6,dtype=None,
                      names=True,skip_header=8)
data7 = np.genfromtxt(filename7,dtype=None,
                      names=True,skip_header=8)
filename2="./IsochroneData/Isochrone2.txt"
filename3="./IsochroneData/Isochrone3.txt"
filename4="./IsochroneData/Isochrone4.txt"
filename5="./IsochroneData/Isochrone5.txt"
filename8="./IsochroneData/Isochrone8.txt"
filename9="./IsochroneData/Isochrone9.txt"
filename12="./IsochroneData/Isochrone12.txt"
filename13="./IsochroneData/Isochrone13.txt"

data2= np.genfromtxt(filename2,dtype=None, names=True,skip_header=8)
data3 = np.genfromtxt(filename3,dtype=None, names=True,skip_header=8)
data4 = np.genfromtxt(filename4,dtype=None, names=True,skip_header=8)
data5 = np.genfromtxt(filename5,dtype=None, names=True,skip_header=8)
data8 = np.genfromtxt(filename8,dtype=None, names=True,skip_header=8)
data9 = np.genfromtxt(filename9,dtype=None, names=True,skip_header=8)
data12 = np.genfromtxt(filename12,dtype=None, names=True,skip_header=8)
data13 = np.genfromtxt(filename13,dtype=None, names=True,skip_header=8)


# In[23]:


# Plot Isochrones 
# For Carina

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot Isochrones

# Isochrone for 1 Gyr
# Plotting Color vs. Difference in Color 
plt.plot(data1['B']-data1['R'], data1['R'], color='blue', 
         linewidth=5, label='1 Gyr')
###EDIT Here, following the same format as the line above 




plt.plot(data2['B']-data2['R'], data2['R'], color='pink', linewidth=5, label='2 Gyr')
plt.plot(data3['B']-data3['R'], data3['R'], color='black', linewidth=5, label='3 Gyr')
plt.plot(data4['B']-data4['R'], data4['R'], color='purple', linewidth=5, label='4 Gyr')
plt.plot(data5['B']-data5['R'], data5['R'], color='orange', linewidth=5, label='5 Gyr')
plt.plot(data6['B']-data6['R'], data6['R'], color='magenta', linewidth=5, label='6 Gyr')
plt.plot(data7['B']-data7['R'], data7['R'], color='green', linewidth=5, label='7 Gyr')
plt.plot(data8['B']-data8['R'], data8['R'], color='brown', linewidth=5, label='8 Gyr')
plt.plot(data9['B']-data9['R'], data9['R'], color='cyan', linewidth=5, label='9 Gyr')
plt.plot(data10['B']-data10['R'], data10['R'], color='red', linewidth=5, label='10 Gyr')
plt.plot(data11['B']-data11['R'], data11['R'], color='yellow', linewidth=5, label='11 Gyr')
plt.plot(data12['B']-data12['R'], data12['R'], color='olive', linewidth=5, label='12 Gyr')
plt.plot(data13['B']-data13['R'], data13['R'], color='gray', linewidth=5, label='13 Gyr')


# Add axis labels
plt.xlabel('B-R', fontsize=22)
plt.ylabel('M$_R$', fontsize=22)

#set axis limits
plt.xlim(-0.5,2)
plt.ylim(5,-2.5)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.figtext(0.5, 0.15, 'CMD for Carina dSph', fontsize=22)
plt.show()
plt.savefig('IsochroneCarina.png')


# # Q2
# 
# Could there be younger ages than suggested in the Tolstoy plot?
# Try adding younger isochrones to the above plot.
# 
# # Q3
# 
# What do you think might cause the bursts of star formation?
# 

# Galactic collisions where nebulae of gas collide to form stars rapidly. Also, supernovae from previous population of stars creating stellar nebulae to increase gas density and decrease temperature.

# In[ ]:




