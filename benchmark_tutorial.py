
# coding: utf-8

# In[1]:

# import general modules
import numpy as np

# import the modules in this package 
import GC_benchmark as GC

# block warnings to make this notebook cleaner.
import warnings
warnings.filterwarnings('ignore')


# In[2]:

# A simple benchmark example. 

# the output files from model1 and model2
filename1='sample_data/GCHP.wetdep.regrid.20130702.nc4' # GCHP output
filename2='sample_data/GEOSCHEM_Diagnostics_Hrly.201307020000.nc' # GC-classic NC diag

# create a "GC_benchmark" object.      
# shortname will be used as the file names of output plots
# longname will be used as the titles of plots. If not specified then use shortname for titles too.
bm = GC.benchmark(outputdir='./outputfig/',shortname='wetdep',longname='Wet Depostion')

# get all tracer data from file1. tag will be used as subplot titles. GCHP output need to be filpped.
# prefix is already set to'SPC_' by default. Just write it out for clarity. Change it when necessary.
# time is already set to 0 by default. Change it to read other time slice.
bm.getdata1(filename1,tag='GCHP',prefix='SPC_',flip=True,time=0)

# now the tracer names and data are hold in this benchmark object.
print(bm.tracername)
print(len(bm.data1))
print(np.shape(bm.data1[0])) # each item in the list is a 3D numpy array (Nlev,Nlat,Nlon)


# In[3]:

# With the knowledge of tracer names, we can read data from file2, matching file1's tracer sequence. 
bm.getdata2(filename2,tag='GCclassic')
print(len(bm.data2))


# In[4]:

# Now the bm object holds all the information needed. A quick benchmark can be done in one line.
bm.plot_all() 

# Some pdfs will be generated in the outputdir specified before. 


# #### They look like that:
# - Just show one page for each file here, but they are multi-page pdfs containing all tracers in bm.tracername
# 
# - For unknown reason, running inside this notebook will stretch the output plot a little bit. Run benchmark_tutorial.py instead. (which is just created by "jupyter nbconvert benchmark_tutorial.ipynb --to python" )
# 
# <img src='outputfig/wetdep_surf.pdf'>
# <img src='outputfig/wetdep_zonalmean.pdf'>

# In[5]:

# delete the benchmark object for the second example.
del bm


# In[6]:

# Now do the benchmark with more steps, to get more control

# those are the same as the previous example
filename1='sample_data/GCHP.wetdep.regrid.20130702.nc4' # GCHP output
filename2='sample_data/GEOSCHEM_Diagnostics_Hrly.201307020000.nc' # GC-classic NC diag
bm = GC.benchmark(outputdir='./outputfig/',shortname='wetdep',longname='Wet Depostion')

# But this time we only read two tracers, NO and O3.
# If tracerlist not specified, then read all.
bm.getdata1(filename1,tag='GCHP',tracerlist=['NO','O3'],
            prefix='SPC_',flip=True,time=0)
# Correspondingly, getdata2 only gets these two tracers, too.
bm.getdata2(filename2,tag='GCclassic')

# now we only have two tracers
print(bm.tracername)
print(len(bm.data1))
print(len(bm.data2))


# In[7]:

# with lower-level methods, you can specify the output filename. 
# Files will still be in the outputdir specified before, but now the shortname doesn't matter

# specify a certain level
bm.plot_layer(pdfname='Only_two_tracer_level8.pdf',lev=7,tag='the 8th level')

# specify a certain range of layers to plot zonal profile 
bm.plot_zonal(pdfname='Only_two_tracer_strato.pdf',mean=True,levs=[32,71],tag='stratosphere zonal mean')

# for more parameters, see help(bm.plot_layer), help(bm.plot_zonal)


# #### The stratospheric one look like that:
# Useful for test StratChem
# <img src='outputfig/Only_two_tracer_strato.pdf'>
