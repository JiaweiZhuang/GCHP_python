"""
PURPOSE:
    GEOSChem/GCHP python benchmark tool. 

NOTES:
    
    1) Assume python3 syntax. Not tested extensively with python2.
    
    2) The key functionality is comparing two models, no matter
       (GCC,GCHP), (GCC,GCC), or (GCHP,GCHP).
       
    3) Currently only benchmark instantaneous tracer mixing ratio.
    
    4) For GEOSChem-classic, only support netCDF files. Will NOT make any
       effort to support bpch format, because NC diagnostics are already  
       implemented in v11-01, and should become standard in v11-02
    
    5) Only use the most basic packages (netCDF4,matplotlib,Basemap) 
       to gain full control over the plot style and format. 

REVISION HISTORY:
    12 Feb 2017 - J.W.Zhuang - Initial version
    
"""

import numpy as np
import matplotlib.pyplot as plt
import gamap as ga
from netCDF4 import Dataset
from matplotlib.backends.backend_pdf import PdfPages
from itertools import chain


class benchmark:
    '''
    Compare two GEOS-Chem NetCDF outputs (either classic or HP), including 
    reading data, making standard plots and quantifying error. 
    '''
    
    def __init__(self,outputdir='./',shortname='DefaultName',longname=None):
        '''
        Just a routinely method to initialize the "benchmark" object.
        Not doing anything meaningful except for setting the case name. 
        [ can be viewed as something like plt.figure() ]
        The most important parameters are further initialized by the 
        getdata1 method.
        
        Parameters
        ----------
        outputdir: string, optional but recommended
            To output directory plots. Use the current dir by default.
            It should end with '/'
            
        shortname: string, optional but recommended
            The shortname of this test. 
            Will be used as the file names of output plots
            
        longname: string, optional but recommended
            The longname of this test. 
            Will be used as the titles of output plots

        '''
        self.outputdir=outputdir
        self.shortname=shortname
        if longname is None:
            self.longname=shortname
        else:
            self.longname=longname
            
        self.isData1=False # haven't read data1
        
    def getdata1(self,filename,tracerlist=None,tag='model1',
                 prefix='SPC_',time=0,flip=False):
        '''
        Get data from the 1st model.
        
        Parameters
        ----------
        filename: string
            The name of the netCDF output file
    
        tracerlist: list of strings, optional
            The names of tracers to be extracted from the netCDF file.
            If not specified, then try to extract all the tracers that 
            matches the prefix.
            
        tag: string, optional but recommended
            The name of the model shown in the plots
            
        prefix: string, optional
            The prefix of variable names in the netCDF file
            
        time: integer, optional
            the time slice to extract. get the first slice by default
        
        flip: logic, optional
            flip the lev dimension or not. Mainly for GCHP I/O issue.
        
        
        Important Returns
        ----------
        self.tracername: list of strings
            tracer names without prefix
        
        self.data1: list of 3D numpy arrays
            3D (lev,lat,lon) data at one time slice.
            
        '''
        
        # open the netcdf file with read-only mode
        fh = Dataset(filename, "r", format="NETCDF4")
        
        # get dimension info
        # use len() instead of .size for compatibility with python2
        self.lon=fh['lon'][:]
        self.lat=fh['lat'][:]
        self.lev=fh['lev'][:]
        self.Nlon=len(self.lon)
        self.Nlat=len(self.lat)
        self.Nlev=len(self.lev)
        
        # initialize an empty list to hold the data.
        self.data1=[]
        
        # use 'is' instead of '==' to test None
        if (tracerlist is None):
            # tracerlist not specified. go through all entries
            
            # initialize an empty list to hold the name
            self.tracername=[]
            
            # get the variable information. fh.variables is a dictionary 
            # consisting of key-value pairs. The key is the variable name 
            # and the value is all the other information
            var = fh.variables
            
            # the standard way to go through a dictionary
            # loop over (key,value) pairs 
            for k,v in var.items():
                # skip the dimension variables
                if prefix in k:
                    # get the name without prefix
                    self.tracername.append(k.replace(prefix,''))
                    # get one time slice
                    data=v[time,:,:,:]

                    if flip: data=data[::-1,:,:]
                    self.data1.append(data)
            
        else:
            # tracerlist specified. 
            self.tracername=tracerlist
            
            for tracer in tracerlist:
                # expand, for example, 'O3' to 'SPC_O3'
                varname=prefix+tracer
                # extract the data directly by the key
                data=fh[varname][time,:,:,:]
    
                if flip: data=data[::-1,:,:]
                self.data1.append(data)

        # always remember to close the NC file
        fh.close() 
        
        # data1 is read
        self.isData1=True
        
        self.model1=tag # for plot names
        
    def getdata2(self,filename,tag='model2',
                 prefix='SPC_',time=0,flip=False):
        '''
        Get data from the 2nd model, with the tracerlist got from getdata1 
        Must run getdata1 first. data2 will match data1's tracer sequence. 
        
        Parameters
        ----------
        see getdata1. The only difference is not requiring tracerlist input
        
        Impoertant Returns
        ----------        
        self.data2: list of 3D numpy arrays
            3D (lev,lat,lon) data at one time slice.
        
        '''
        if self.isData1 == False:
            raise ValueError('must run getdata1 first')
            
        # open the netcdf file with read-only mode
        fh = Dataset(filename, "r", format="NETCDF4")
            
        # check dimension 
        if self.Nlon != len(fh['lon'][:]):
            raise ValueError('lon dimension does not match')
        if self.Nlat != len(fh['lat'][:]):
            raise ValueError('lat dimension does not match')
        if self.Nlev != len(fh['lev'][:]):
            raise ValueError('lev dimension does not match')

        # initialize an empty list to hold the data.
        self.data2=[]
            
        for tracer in self.tracername:
            # expand, for example, 'O3' to 'SPC_O3'
            varname=prefix+tracer
            # extract the data directly by the key
            data=fh[varname][time,:,:,:]
    
            if flip: data=data[::-1,:,:]
            self.data2.append(data)
            
        # always remember to close the NC file
        fh.close() 
        
        self.model2=tag # for plot names
        
    def getdata0(self,filename,
                 prefix='SPC_',time=0,flip=False):
        '''
        Get initial condition, with the tracerlist got from getdata1 
        Must run getdata1 first. data0 will match data1's tracer sequence. 
        
        Parameters
        ----------
        see getdata1. The only difference is not requiring tracerlist input
        
        Impoertant Returns
        ----------        
        self.data0: list of 3D numpy arrays
            3D (lev,lat,lon) data at one time slice.
        
        '''
        if self.isData1 == False:
            raise ValueError('must run getdata1 first')
            
        # open the netcdf file with read-only mode
        fh = Dataset(filename, "r", format="NETCDF4")
            
        # check dimension 
        if self.Nlon != len(fh['lon'][:]):
            raise ValueError('lon dimension does not match')
        if self.Nlat != len(fh['lat'][:]):
            raise ValueError('lat dimension does not match')
        if self.Nlev != len(fh['lev'][:]):
            raise ValueError('lev dimension does not match')

        # initialize an empty list to hold the data.
        self.data0=[]
            
        for tracer in self.tracername:
            # expand, for example, 'O3' to 'SPC_O3'
            varname=prefix+tracer
            # extract the data directly by the key
            data=fh[varname][time,:,:,:]
    
            if flip: data=data[::-1,:,:]
            self.data0.append(data)
            
        # always remember to close the NC file
        fh.close() 
      
    def plot_layer(self,lev=0,plot_change=False,switch_scale=False,
                   pdfname='default_layerplot.pdf',tag='',rows=3):
        '''
        Compare self.data1 and self.data2 on a specific level. Loop over all
        tracers in self.tracerlist. Plot on one multi-page pdf file. 
        
        Parameters
        ----------
        lev: integer
            The level to compare
            
        pdfname: string, optional but recommended
            The name the output pdf file
            
        tag: string, optional but recommended
            Will be shown as part of the title. For example,'surface', '500hpa'
            
        rows: integer, optional
            How many rows on a page. Although the page size will be 
            adjusted accordingly, 3 rows look better.
            
        plot_change: logical, optional
            If set to True, then plot the change with respect to the initial
            condition (i.e. self.data0), rather than the original field.
            In this case, self.getdata0 must be executed before.
        
        switch_scale: logical, optional
            By default, the scale of data2 is used for both data1 and data2.
            This is because model2 is often regarded as the more correct one
            (the one for reference), so if model1 goes crazy, the color scale 
            will still make sense.
            
            If set to True, then use the scale of data1 instead. This is often
            combined with plot_change=True to more clearly show the regridding
            error.
            
    
        Important Returns
        ----------
            A pdf file containing the plots of all tracers on that level.
        
        '''
        
        N = len(self.tracername) # number of tracers
        Npages = (N-1)//rows+1 # have many pdf pages needed
        
        print('create ',pdfname)
        pdf=PdfPages(self.outputdir+pdfname) # create the pdf to save figures
        
        for ipage in range(Npages):
            
            fig, axarr = plt.subplots(rows,3,figsize=(12,rows*3))
        
            for i in range(rows):
                i_tracer = ipage*rows + i
                print('plotting: no.',i_tracer+1,self.tracername[i_tracer])
                
                # make variable names shorter for simplicity
                tracername = self.tracername[i_tracer]
                
                # assume v/v, convert to ppbv
                # might need modification in the future
                data1 = self.data1[i_tracer][lev,:,:]*1e9
                data2 = self.data2[i_tracer][lev,:,:]*1e9
                unit='ppbv'
                
                if plot_change:
                    # plot the change with respect to the initial condition,
                    # instead to the original field
                    data0 = self.data0[i_tracer][lev,:,:]*1e9
                    data1 -= data0
                    data2 -= data0
                    cmap=plt.cm.RdBu_r
                    
                else:
                    cmap=ga.WhGrYlRd
                
                if np.max(np.abs(data1)) < 1e-1 :
                    data1 *= 1e3
                    data2 *= 1e3
                    unit='pptv'             
                elif np.max(np.abs(data1)) > 1e3 :
                    data1 /= 1e3
                    data2 /= 1e3
                    unit='ppmv'
                                  
                # use the same scale for data1 and data2
                # use the scale of data2 by default.
                if switch_scale:
                    range_data = np.max(np.abs(data1))
                else:
                    range_data = np.max(np.abs(data2))
                # calculate the difference between two data sets
                data_diff = data1-data2
                range_diff=np.max(np.abs(data_diff))
                
                if plot_change:
                    vmin = -range_data
                    vmax = range_data
                    title= tracername+' change; '
                else:
                    vmin = 0
                    vmax = range_data
                    title= tracername+'; '
                            
                ga.tvmap(data1,axis=axarr[i,0],vmin=vmin,vmax=vmax,unit=unit,
                         cmap=cmap,title=title+self.model1,ticks = False)
                ga.tvmap(data2,axis=axarr[i,1],vmin=vmin,vmax=vmax,unit=unit,
                         cmap=cmap,title=title+self.model2,ticks = False)
                ga.tvmap(data_diff,axis=axarr[i,2],unit=unit,
                         title=self.model1+' — '+self.model2,ticks = False,
                         cmap=plt.cm.RdBu_r,vmax=range_diff,vmin=-range_diff)
                
                if i_tracer+1 == N : 
                    i_hide = rows-1
                    while(i_hide > i):
                        [a.axis('off') for a in axarr[i_hide, :]]
                        i_hide -= 1
                    break # not using the full page
                
            fig.suptitle(self.longname+'; '+tag,fontsize=15)

            print('saving one pdf page')
            pdf.savefig(fig)  # saves the figure into a pdf page
            plt.close() # close the current figure, prepare for the next page
            
        pdf.close() # close the pdf after saving all figures
        print(pdfname,' finished')
        
    def plot_zonal(self,mean=False,ilon=0,levs=None,
                   plot_change=False,switch_scale=False,
                   pdfname='default_zonalplot.pdf',tag='',rows=3):
        
        '''
        Compare the zonal profiles of self.data1 and self.data2. Loop over all
        tracers in self.tracerlist. Plot on one multi-page pdf file. 
        
        Parameters
        ----------
        mean: logical, optional
            If set to True, then plot the zonal mean.
            Otherwise plot one cross-section specified by ilon
            
        ilon: integer, optional
            The longitude index of the cross-section. 
            Will have not effect mean is True.
            
        pdfname: string, optional but recommended
            The name the output pdf file
            
        tag: string, optional but recommended
            Will be shown as part of the title. For example,'zonal mean'
            
        rows: integer, optional
            How many rows on a page. Although the page size will be 
            adjusted accordingly, 3 rows look better.
            
        plot_change: logical, optional
            If set to True, then plot the change with respect to the initial
            condition (i.e. self.data0), rather than the original field.
            In this case, self.getdata0 must be executed before.
        
        switch_scale: logical, optional
            By default, the scale of data2 is used for both data1 and data2.
            This is because model2 is often regarded as the more correct one
            (the one for reference), so if model1 goes crazy, the color scale 
            will still make sense.
            
            If set to True, then use the scale of data1 instead. This is often
            combined with plot_change=True to more clearly show the regridding
            error.
    
        Important Returns
        ----------
            A pdf file containing the plots of all tracers' zonal profile.
        
        '''
                
        N = len(self.tracername) # number of tracers
        
        Npages = (N-1)//rows+1 # have many pdf pages needed
        
        print('create ',pdfname)
        pdf=PdfPages(self.outputdir+pdfname) # create the pdf to save figures
        
        for ipage in range(Npages):
            
            fig, axarr = plt.subplots(rows,3,figsize=(12,rows*3))
        
            for i in range(rows):
                i_tracer = ipage*rows + i
                print('plotting: no.',i_tracer+1,self.tracername[i_tracer])
                
                # make variable names shorter for simplicity.
                # '*1.0' is needed, otherwise data1_3D amd self.data1 will 
                # share the same physical address, and will affect next plots.
                tracername = self.tracername[i_tracer]
                data1_3D = self.data1[i_tracer][:]*1.0
                data2_3D = self.data2[i_tracer][:]*1.0
                                     
                if plot_change:
                    # plot the change with respect to the initial condition,
                    # instead to the original field
                    data0_3D = self.data0[i_tracer][:]*1.0
                    data1_3D -= data0_3D
                    data2_3D -= data0_3D
                    cmap=plt.cm.RdBu_r
                else:
                    cmap=ga.WhGrYlRd
                
                # get limited levels if requested
                if levs is None:
                    lev = self.lev
                else:
                    lev = self.lev[levs[0]:levs[1]]
                    data1_3D=data1_3D[levs[0]:levs[1],:,:]
                    data2_3D=data2_3D[levs[0]:levs[1],:,:]
                
                if mean:
                    # get zonal mean
                    data1 = np.mean(data1_3D,axis=2)
                    data2 = np.mean(data2_3D,axis=2)
                else:
                    # get cross-section
                    data1 = data1_3D[:,:,ilon]
                    data2 = data2_3D[:,:,ilon]
                    
                # assume v/v, convert to ppbv
                # might need modification in the future
                data1 *= 1e9
                data2 *= 1e9
                unit='ppbv'

                if np.max(data1) < 1e-1 :
                    data1 *= 1e3
                    data2 *= 1e3
                    unit='pptv'                    
                elif np.max(data1) > 1e3 :
                    data1 /= 1e3
                    data2 /= 1e3
                    unit='ppmv'
                    
                # use the same scale for data1 and data2
                # use the scale of data2 by default.
                if switch_scale:
                    range_data = np.max(np.abs(data1))
                else:
                    range_data = np.max(np.abs(data2))
                # calculate the difference between two data sets
                data_diff = data1-data2
                range_diff=np.max(np.abs(data_diff))
                
                if plot_change:
                    vmin = -range_data
                    vmax = range_data
                    title= tracername+' change; '
                else:
                    vmin = 0
                    vmax = range_data
                    title= tracername+'; '
                         
                xlabel='lat'
                ylabel='level'
                            
                ga.tvplot(data1,axis=axarr[i,0],vmin=vmin,vmax=vmax,unit=unit,
                          cmap=cmap,x=self.lat,y=lev,xlabel=xlabel,ylabel=ylabel,
                          title=title+self.model1)
                ga.tvplot(data2,axis=axarr[i,1],vmin=vmin,vmax=vmax,unit=unit,
                          cmap=cmap,x=self.lat,y=lev,xlabel=xlabel,ylabel=ylabel,
                          title=title+self.model2)
                ga.tvplot(data_diff,axis=axarr[i,2],unit=unit,
                          x=self.lat,y=lev,xlabel=xlabel,ylabel=ylabel,
                          title=self.model1+' — '+self.model2,
                          cmap=plt.cm.RdBu_r,vmax=range_diff,vmin=-range_diff)
                
                # hide x ticks except the bottom plots
                if (i != rows-1 ) and ( i_tracer+1 != N) :
                    plt.setp([a.get_xticklabels() for a in axarr[i, :]], 
                             visible=False)
                    plt.setp([a.set_xlabel('') for a in axarr[i, :]], 
                             visible=False)
                
                if i_tracer+1 == N : 
                    if i != rows-1 :
                        [a.axis('off') for a in axarr[i+1, :]]
                    break # not using the full page
                
            #  hide y ticks for right plots
            plt.setp([a.get_yticklabels() for a in list(chain.from_iterable(axarr[:, 1:3]))], 
                     visible=False)
            plt.setp([a.set_ylabel('') for a in list(chain.from_iterable(axarr[:, 1:3]))],
                     visible=False)
                
            fig.suptitle(self.longname+'; '+tag,fontsize=15)

            print('saving one pdf page')
            pdf.savefig(fig)  # saves the figure into a pdf page
            plt.close() # close the current figure, prepare for the next page
            
        pdf.close() # close the pdf after saving all figures
        print(pdfname,' finished')
        
        
    def plot_all(self,plot_surf=True,plot_500hpa=True,plot_zonal=True,
                 plot_change=False):
        '''
        Just to wrap several plot_layer and plot_zonal calls for convenience.
        
        Parameters
        ----------
        plot_change: logical, optional
            If set to True, then plot the change with respect to the initial
            condition (i.e. self.data0), rather than the original field.
            In this case, self.getdata0 must be executed before.
        
        plot_surf,plot_500hpa,plot_zonal: logical, optional
            plot that feature or not. Default is to plot everything.
        
        Important Returns
        ----------
            Several pdf files.
        '''
        
        if plot_change:
            switch_scale=True # use data1's scale when plotting change
            tag='_change' # the file name should be slightly different
        else:
            switch_scale=False
            tag=''
        
        if plot_surf:
            self.plot_layer(pdfname=self.shortname+tag+'_surf.pdf',lev=0,tag='surface',
                            switch_scale=switch_scale,plot_change=plot_change)
        if plot_500hpa:
            self.plot_layer(pdfname=self.shortname+tag+'_500hpa.pdf',lev=22,tag='500hpa',
                            switch_scale=switch_scale,plot_change=plot_change)
        if plot_zonal:
            self.plot_zonal(pdfname=self.shortname+tag+'_180lon.pdf',
                            ilon=0,levs=[0,20],tag='180lon',
                            switch_scale=switch_scale,plot_change=plot_change)
            self.plot_zonal(pdfname=self.shortname+tag+'_zonalmean.pdf',
                            mean=True,levs=[0,20],tag='zonal mean',
                            switch_scale=switch_scale,plot_change=plot_change)
            

                
