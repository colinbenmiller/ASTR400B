#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import astropy.units as u


# In[2]:


def Read(filename): 
    """
    this function reads in the filename
    input is the filename
    output is the time of the dataset, the set of particles, and data associated with it
    """
    file = open(filename, 'r' ) #open file
    line1 = file.readline() #read first line of file
    label, value = line1.split() #splits the lines of up to retrieve
    time = float(value)*u.Myr #intializes the time variable
    line2 = file.readline() #reads the second line
    label, value = line2.split() #splits the lines up
    particles = int(value) #intializes the particle variable output
    file.close() #closes the file
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3) #keeps the column header amd starts on line 4 for the data and generates data
    return time, particles, data #returns outputs


# In[3]:


filename = "./MW_000.txt" #this brings and initializes the file


# In[4]:


Read(filename) #this reads the file


# In[12]:


print(data['type'][1]) #this tests the file


# In[ ]:




