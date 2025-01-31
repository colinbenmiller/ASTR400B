#!/usr/bin/env python
# coding: utf-8

# In[14]:


filename = './MW_000.txt' #initalizes the file


# In[15]:


from ReadFile import Read #imports files
from ParticleProperties import ParticleInfo
import astropy.units as u


# In[16]:


#1 is dark matter, 2 is disk, 3 is bulge


# In[17]:


ParticleInfo(filename, 2, 100) #prints out specific data


# In[18]:


threeDdistance, threeDvelocity, mass = ParticleInfo(filename, 2, 100) #separates out the variables with new names


# In[19]:


threeDdistance.to(u.lightyear) #converts to lightyears


# In[ ]:




