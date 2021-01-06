#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Necessary libraries
#from IPython.display import clear_output

#from tqdm import tqdm_notebook as tqdm
import synapseclient
import pandas as pd
import subprocess
import os.path
import getpass

# Set pathway
local = "/data/home/diegocoelho/SynapseConsortium/ROSMap"


# In[4]:


# Import dataframe containing all .fastq files
tablePath = local + "/data/ROSMapRNAseq_DLPFC_table.csv"
df = pd.read_csv(tablePath)
df.head() # visualize


# In[5]:


# Connect on server
syn = synapseclient.login('diegomscoelho','13944066') # Please do not share user content


# In[6]:


# Download and Realignment

df = df.sort_values(by=["name"]) # Sort dataframe by name
samples = df.specimenID.unique() # Samples list

fqPath = local + "/fastq"
klPath = local + "/kallisto/"
os.popen('mkdir '+fqPath) # Create fastq Path
os.popen('mkdir '+klPath) # Create kallisto output path

# Directory to builded kallisto index
kalRef = "/data/home/diegocoelho/refs/human_kallisto/Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx"

def runKallisto(samples,path):
    
    # Create sample directory path
    kalDir = klPath + samples
    os.popen('mkdir '+kalDir)
    
    # Import ids
    
    print("\nSamples:")
    name1, name2 = df.name[df.specimenID==samples]
    r1, r2 = df.id[df.specimenID==samples]
    print("Sample.1: "+name1)
    print("r1: "+r1)
    print("Sample.2: "+name2)
    print("r2: "+r2)
    
    if os.path.isfile(kalDir+'/abundance.tsv'):
        print("Sample " + samples + " was already processed!")
        pass
    
    else:
        print("Downloading " + samples + " ...")
        file1 = syn.get(r1, downloadLocation=fqPath) # Download file1 to kallisto
        file2 = syn.get(r2, downloadLocation=fqPath) # Download file2 to kallisto

        # Run Kallisto
        print("Running kallisto on " + samples + " ...")
        subprocess.call('kallisto quant -i '+
                 kalRef + ' -t 7 -l 100 -s 20 -o '+ kalDir+ ' '+
                file1.path+' '+file2.path, shell=True)

        # If necessary, remove files after using
        os.popen('rm '+file1.path)
        os.popen('rm '+file2.path)


# In[ ]:


#get_ipython().run_line_magic('timeit', '')
# Run process
for i in samples:
    runKallisto(samples = i, path = local)


# In[ ]:




