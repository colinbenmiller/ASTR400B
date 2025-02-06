#!/usr/bin/env python
# coding: utf-8

# In[4]:


from ReadFile import Read #imports read file, numpy, and astropy units
import numpy as np
import astropy.units as u


# In[7]:


def ComponentMass(filename, particletype):
    """ 
    this function adds reads in the filename and particle type and intitialzes the variables to be used to then find the specific data
    in the file and calculates the component section of galaxy's mass and in terms of mass of the sun
    """
    time, particles, data = Read(filename) #separates file data into arrays to store data to specific variables 
    index=np.where(data['type']==particletype) #creates an index to find specific line of data and renames it to input parameter    
    masses = data['m'][index] #intializes value of mass to data from mass column based on particlerow input
    massT = 0 #initializes to zero
    massT = np.sum(masses) #sum array for particular particle type
    massT = massT* (10**(-2))*u.Msun #makes in unit of solar masses
    massT = np.round(massT, 3) #rounds mass total to 3 decimals
    return massT #returns calculated values
    


# In[ ]:




