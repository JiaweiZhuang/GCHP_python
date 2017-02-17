#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PURPOSE:
    gamap-like routines for plotting, to ease the transition from 
    IDL/gamap to python.
    
NOTES:
    1) The WhGrYlRd scheme data is taken from
       https://github.com/thackray/geosplot/issues/4
    
REVISION HISTORY:
    12 Feb 2017 - J.W.Zhuang - Initial version
    
"""
import matplotlib.pyplot as plt
import numpy as np
import os
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap

# get gamap's WhGrYlRd color scheme from file
current_dir = os.path.dirname(os.path.abspath(__file__))
WhGrYlRd_scheme = np.genfromtxt(current_dir+'/WhGrYlRd.txt',delimiter=' ')
WhGrYlRd = ListedColormap(WhGrYlRd_scheme/255.0)

def tvmap(data,axis=None,title='',
          lon=None,lat=None,vmax=None,vmin=None,unit='unit',
          continent = True, grid = True, ticks=True,
          cmap=WhGrYlRd):
    '''
    Plot a (geographical) 2D array on a map, very similar to 
    IDL/gamap's tvmap routine.
    
    Parameters
    ----------
    data: 2D numpy array
        The data to plot

    axis: matplotlib axis object, optional
        The axis to plot on. If not specified, then create a new plot and
        shows the image immediately.
        
    title: string, optional
        title shown on the top of the plot
        
    lon: 1D numpy array, optional
        The longitude range of the data. Assume 180W~180E by default
        
    lat: 1D numpy array, optional
        The latitude range of the data. 
        Assume 90S~90N with half-polar cell by default
        
    cmap: matplotlib colormap object, optional
        The colormap to use. By default, use gamap's 'WhGrYlRd' scheme, 
        good for showing tracer concentration. 
        To plot difference, plt.cm.RdBu_r is recommended.
    
    vmax, vmin: real, optional
        The max/min of colorbar. Use the range of data by default.
    
    unit: string, optional
        unit shown near the colorbar
        
    continent: logical, optional
        show the continent boundary or not
    
    grid, ticks:  logical, optional
        show the grid line (and ticks) or not

    Important Returns
    ----------
        A plot on the map
        
    '''
    
    # check dimension first
    if np.ndim(data) != 2:
        raise ValueError('input data should be a 2D array!')
        
    # the shape of data should be (Nlat,Nlon)
    Nlat, Nlon = np.shape(data)
    
    # assume a global grid with half-polar cell if lon and lat are not specified
    if lon is None:
        lon=np.linspace(-180,180,Nlon,endpoint=False)
    if lat is None:
        lat=np.linspace(-90,90,Nlat,endpoint=True)

    # center to bound, needed by pcolormesh
    lon_b = 0.5*(lon[1:]+lon[:-1]) # inner bound
    lon_b = np.append(2*lon_b[0]-lon_b[1],lon_b) # left bound
    lon_b = np.append(lon_b,2*lon_b[-1]-lon_b[-2]) # right bound
    
    lat_b = 0.5*(lat[1:]+lat[:-1]) # inner bound
    lat_b = np.append( max(-90,2*lat_b[0]-lat_b[1]) , lat_b  ) # left bound
    lat_b = np.append( lat_b , min(90,2*lat_b[-1]-lat_b[-2]) ) # right bound
    
    # 2D grid mesh
    lons, lats = np.meshgrid(lon_b,lat_b)
    
    # If axis is not specified, create a new one. Otherwise plot on the given axis
    show = False
    if axis is None: 
        plt.figure(figsize=(12,6))
        axis = plt.gca()
        # If axis is not specified, show the plot immediately 
        # i.e. "tvmap in one line"
        show = True

    # 'cyl' projection plot all boxes at the same size, as in IDL/gamap
    m = Basemap(ax=axis, projection='cyl',
                llcrnrlat=lat_b[0],urcrnrlat=lat_b[-1],
                llcrnrlon=lon_b[0],urcrnrlon=lon_b[-1],
                resolution='c',fix_aspect=False)
    
    if continent : m.drawcoastlines()
    
    if grid :
        # plot about 6 grid lines. Also work for regional grid.
        lon_stride = 60/np.round(360/(lon[-1]-lon[0]))
        lat_stride = 30/np.round(180/(lat[-1]-lat[0]))
        if ticks: 
            labels=[1,0,0,0]
        else: 
            labels=[0,0,0,0]
        m.drawparallels(np.arange(-90,90,lat_stride),labels=labels)
        if ticks: 
            labels=[0,0,0,1]
        else: 
            labels=[0,0,0,0]
        m.drawmeridians(np.arange(-180,180,lon_stride),labels=labels)

    # linewidth=0 and rasterized=True remove white lines between boxes, 
    # especially for pdf and postscript formats
    im = m.pcolormesh(lons,lats,data,cmap=cmap,vmax=vmax,vmin=vmin,
                      linewidth=0,rasterized=True)
    cb = m.colorbar(im)
    cb.ax.set_title(unit,fontsize=12,y=1.0,x=2.0) # on the top of the color bar
    plt.title(title)
    
    if show: plt.show()
    
def tvplot(data,axis=None,title='',xlabel='x',ylabel='y',
           x=None,y=None,vmax=None,vmin=None,unit='unit',
           cmap=WhGrYlRd):
    '''
    Plot a general 2D array, very similar to IDL/gamap's tvplot routine.
    
    Parameters
    ----------
    data: 2D numpy array
        The data to plot

    axis: matplotlib axis object, optional
        The axis to plot on. If not specified, then create a new plot and
        shows the image immediately.
        
    title: string, optional
        title shown on the top of the plot
        
    x: 1D numpy array, optional
        x range of the data. Use array index by default
        
    y: 1D numpy array, optional
        y range of the data. Use array index by default
        
    cmap: matplotlib colormap object, optional
        The colormap to use. By default, use gamap's 'WhGrYlRd' scheme, 
        good for showing tracer concentration. 
        To plot difference, plt.cm.RdBu_r is recommended.
    
    vmax, vmin: real, optional
        The max/min of colorbar. Use the range of data by default.
    
    unit: string, optional
        unit shown near the colorbar
        
    xlabel,ylabel: string, optional
        labels for x and y dimensions
        
    Important Returns
    ----------
        A general 2D plot without a the map
        
    '''

    # the shape of data should be (Ny,Nx)
    Ny, Nx = np.shape(data)
    
    # assume integer labeling if x and y are not specified
    if x is None:
        x=np.arange(1,Nx+1)
    if y is None:
        y=np.arange(1,Ny+1)

    # center to bound, needed by pcolormesh
    x_b = 0.5*(x[1:]+x[:-1]) # inner bound
    x_b = np.append(2*x_b[0]-x_b[1],x_b) # left bound
    x_b = np.append(x_b,2*x_b[-1]-x_b[-2]) # right bound
    
    y_b = 0.5*(y[1:]+y[:-1]) # inner bound
    y_b = np.append( 2*y_b[0]-y_b[1] ,  y_b  ) # left bound
    y_b = np.append( y_b , 2*y_b[-1]-y_b[-2] ) # right bound
    
    # 2D grid mesh
    xs, ys = np.meshgrid(x_b,y_b)

    # If axis is not specified, create a new one. Otherwise plot on the given axis
    show = False
    if axis is None: 
        plt.figure(figsize=(12,6))
        axis = plt.gca()
        # If axis is not specified, show the plot immediately 
        # i.e. "tvmap in one line"
        show = True
        
    im = axis.pcolormesh(xs,ys,data,cmap=cmap,vmin=vmin,vmax=vmax,
                         linewidth=0,rasterized=True)
    axis.set_xlim( [min(x_b),max(x_b)] )
    axis.set_ylim( [min(y_b),max(y_b)] )
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    cb = plt.colorbar(im,ax=axis)
    cb.ax.set_title(unit,fontsize=12,y=1.0,x=2.0) # on the top of the color bar
    axis.set_title(title)
    
    if show: plt.show()
    
    
