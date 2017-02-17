
# coding: utf-8

# #### See benchmark_tutorial.py for the pure python code (which is just created by "jupyter nbconvert benchmark_tutorial.ipynb --to python")

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
# Typically, model2 is regarded as a reference (e.g. old version), 
# and model1 is what you want to examine.
# the sequence will affect the choice of the color scale.
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


# #### View the output pdfs here
# 
# [wetdep_surf.pdf](outputfig/wetdep_surf.pdf)
# 
# [wetdep_500hpa.pdf](outputfig/wetdep_500hpa.pdf)
# 
# [wetdep_180lon.pdf](outputfig/wetdep_180lon.pdf)
# 
# [wetdep_zonalmean.pdf](outputfig/wetdep_zonalmean.pdf)

# In[5]:

# For a single component test, if makes more sense to 
# plot the change in the tracer field instead of the tracer field itself.

# here we extract 5 tracers from the standard restart file to make the size small
# (by "ncks -v SPC_NO,SPC_O3,SPC_NIT,SPC_HNO3,SPC_CH2O gcc_restart.nc gcc_restart_5tracer.nc")
# Indeed it also works with the standard restart file
filename0='sample_data/gcc_restart_5tracer.nc'

# you don't need to run getdata0 if you don't need to plot the change
bm.getdata0(filename0,prefix='SPC_',flip=False)
print(len(bm.data0)) # now bm holds the initial condition


# In[6]:

# set plot_change to True to plot the change to initial condition.
# make sure data0 has already been read in.
bm.plot_all(plot_change=True) 


# #### View the output pdfs here
# 
# [wetdep_change_surf.pdf](outputfig/wetdep_change_surf.pdf)
# 
# [wetdep_change_500hpa.pdf](outputfig/wetdep_change_500hpa.pdf)
# 
# [wetdep_change_180lon.pdf](outputfig/wetdep_change_180lon.pdf)
# 
# [wetdep_change_zonalmean.pdf](outputfig/wetdep_change_zonalmean.pdf)
# 
# - Here we use GCC's 4x5 initial condition, so the "change" in O3 and NO is just the regridding error, because they are not removed by wet deposition. To make a fair comparision, we might use a finer grid as input and final output.

# In[7]:

# delete the benchmark object for the second example.
del bm


# In[8]:

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


# In[9]:

# with lower-level methods, you can specify the output filename. 
# Files will still be in the outputdir specified before, but now the shortname doesn't matter

# specify a certain level
bm.plot_layer(pdfname='Only_two_tracer_level8.pdf',lev=7,tag='the 8th level')

# specify a certain range of layers to plot zonal profile 
bm.plot_zonal(pdfname='Only_two_tracer_strato.pdf',mean=True,levs=[32,71],tag='stratosphere zonal mean')

# plot the change if you've read the initial condition.
filename0='sample_data/gcc_restart_5tracer.nc'
bm.getdata0(filename0)
bm.plot_zonal(pdfname='Only_two_tracer_change.pdf',plot_change=True,switch_scale=True,
              ilon=0,levs=[0,15],tag='lower troposphere 180lon')
'''
note: 
when plotting the change, it is recommended to set switch_scale to True, 
to more clearly show the regridding error in GCHP. 
It ensures GCHP's data range is used for the common color bar.

when plotting the original field, GC-classic(model2)'s color scale is used,
so the plot will still make sense it model1 goes crazy.
'''
# for more parameters, see help(bm.plot_layer), help(bm.plot_zonal)


# #### View the output pdfs here
# 
# [Only_two_tracer_level8.pdf](outputfig/Only_two_tracer_level8.pdf)
# 
# [Only_two_tracer_strato.pdf](outputfig/Only_two_tracer_strato.pdf)
# 
# [Only_two_tracer_change.pdf](outputfig/Only_two_tracer_change.pdf)
