#!/usr/bin/env python
# coding: utf-8

# In[11]:


from ReadFile import Read #imports read file, numpy, and astropy units
import numpy as np
import astropy.units as u


# In[10]:


def ParticleInfo(filename, particletype, particlenumber):
    """ 
    this function adds reads in the filename, particle type, and particle number and intitialzes the variables to be used to then find the specific data
    in the file and calculates 3D velocity in km/s, 3D distance from galactic center in kpc, and mass and in terms of mass of the sun
    """
    time, particles, data = Read(filename) #separates file data into arrays to store data to specific variables 
    index=np.where(data['type'])==particletype #creates an index to find specific line of data and renames it to input parameter
    typenew = data[index] #renames data array with specific line in the file and stores it
    particlerow = typenew[particlenumber -1] #renames value from previous to line up with input paramter and subtracts 1 due to first line being 0
    x=particlerow['x'] #intializes x distance coord position in kpc from said specific column based on particle row
    y=particlerow['y'] #intializes y distance coord position in kpc from said specific column based on particle row
    z=particlerow['z'] #intializes z distance coord position in kpc from said specific column based on particle row
    d=np.sqrt(x**2 + y**2 + z**2) #creates new function to find the magnitude of the 3D position data based on inputs
    d=np.round(d,3) #rounds value to 3 decimal points
    d=d*u.kpc #turns data in kpc to read out when running code
    vx=particlerow['vx'] #intializes x velocity coord position in km/s from said specific column based on particle row
    vy=particlerow['vy'] #intializes y velocity coord position in km/s from said specific column based on particle row
    vz=particlerow['vz'] #intializes z velocity coord position in km/s from said specific column based on particle row
    v=np.sqrt(vx**2 + vy**2 + vz**2) #adds up 3D velocity magnitude based on inputs of particle row data
    v=np.round(v,3) #rounds value to 3 decimal points
    v = v*(u.km/u.s) #turns data in km/s to read out when running code
    mass = particlerow['m'] #intializes value of mass to data from mass column based on particlerow input
    mass = mass * (10**10)*u.Msun #converts mass given from input in terms of mass of the sun
    return d, v, mass #returns calculated values


# In[ ]:




