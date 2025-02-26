#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Homework 6
# Colin Miller


# In[2]:


# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass



# In[3]:


def OrbitCom(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs: takes in galaxy name, snapshot to start, snapshot to end on, and n integer intervals returning COM
    outputs: time and COM position and velocity vectors of a given galaxy in each snapshot 
    """
    
    # compose the filename for output
    fileout = f'Orbit_{galaxy}.txt'
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1
    if galaxy == "M33":
        volDec = 4
    else:
        volDec = 2
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end+1, n)
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids), 7])
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        ilbl = '000' + str(snap_id)
        ilbl = ilbl[-3:]
        filename = f'{galaxy}_{ilbl}.txt'
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename, 2)
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        com_pos = COM.COM_P(delta = delta, volDec = volDec)
        com_vel = COM.COM_V(com_pos[0], com_pos[1], com_pos[2])
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        orbit[i,0] = COM.time.to(u.Gyr).value
        orbit[i,1] = com_pos[0].value
        orbit[i,2] = com_pos[1].value
        orbit[i,3] = com_pos[2].value
        orbit[i,4] = com_vel[0].value
        orbit[i,5] = com_vel[1].value
        orbit[i,6] = com_vel[2].value
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
    print(f"Orbit saved in {fileout}")


# In[4]:


# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
MW = OrbitCom("MW", 0, 800, 5)
M31 = OrbitCom("M31", 0, 800, 5)
M33 = OrbitCom("M33", 0, 800, 5)


# In[5]:


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
data_MW = np.genfromtxt("Orbit_MW.txt", comments = '#')
data_M31 = np.genfromtxt("Orbit_M31.txt", comments = '#')
data_M33 = np.genfromtxt("Orbit_M33.txt", comments = '#')


# In[6]:


# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
tMW = data_MW[:, 0]
xMW = data_MW[:, 1]
yMW = data_MW[:, 2]
zMW = data_MW[:, 3]
vxMW = data_MW[:, 4]
vyMW = data_MW[:, 5]
vzMW = data_MW[:, 6]

tM31 = data_M31[:, 0]
xM31 = data_M31[:, 1]
yM31 = data_M31[:, 2]
zM31 = data_M31[:, 3]
vxM31 = data_M31[:, 4]
vyM31 = data_M31[:, 5]
vzM31 = data_M31[:, 6]

tM33 = data_M33[:, 0]
xM33 = data_M33[:, 1]
yM33 = data_M33[:, 2]
zM33 = data_M33[:, 3]
vxM33 = data_M33[:, 4]
vyM33 = data_M33[:, 5]
vzM33 = data_M33[:, 6]


# In[7]:


# Determine the magnitude of the relative position and velocities 
def MagVectorDiff(x, y, z, x1, y1, z1):
    "inputs: xyz coordinates of each galaxy used"
    "outputs: magnitude of computed vector differnece"
    vmag = np.sqrt((x1 - x)**2 + (y1 - y)**2 + (z1 - z)**2)
    return vmag
# of MW and M31
p_MW_M31 = MagVectorDiff(xMW, yMW, zMW, xM31, yM31, zM31)
v_MW_M31 = MagVectorDiff(vxMW, vyMW, vzMW, vxM31, vyM31, vzM31)                      
# of M33 and M31
p_M33_M31 = MagVectorDiff(xM33, yM33, zM33, xM31, yM31, zM31)
v_M33_M31 = MagVectorDiff(vxM33, vyM33, vzM33, vxM31, vyM31, vzM31)   


# In[19]:


fig, ax = plt.subplots(1, 2, figsize = (25,10))
# Plot the Orbit of the galaxies 
#################################
ax[0].plot(tMW, p_MW_M31, color = 'brown', label = 'MW-M31')
ax[0].plot(tM31, p_M33_M31, color = 'green', label = 'M33-M31')
ax[0].set_xlabel('Time (Gyr)')
ax[0].set_ylabel('Separation (kpc)')
ax[0].set_title('Separation vs. Time')
ax[0].legend()

# Plot the orbital velocities of the galaxies 
#################################
ax[1].plot(tMW, v_MW_M31, color = 'blue', label = 'MW-M31')
ax[1].plot(tM31, v_M33_M31, color = 'black', label = 'M33-M31')
ax[1].set_xlabel('Time (Gyr)')
ax[1].set_ylabel('Rel Velocity (km/s)')
ax[1].set_title('Rel Veloicty vs. Time')
ax[1].legend()

plt.tight_layout()
plt.show()


# In[ ]:





# 1. They will have 2 close encounters before colliding.

# 2. When separation peaks(i.e. the first derivative is zero), the relative velocity hits a low point where the first derivative is zero. This relationship holds true for the inverse when separation hits lowest on curve, the relative velocity peaks.

# 3. MW and M31 merge around 6.75 Gyrs. When this happens, M33's orbital velocity gets smaller and smaller on the peaks and troughs around M31 and gets closer to M31 over time, orbiting elliptically.

# In[ ]:




