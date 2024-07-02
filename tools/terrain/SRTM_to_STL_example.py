def SRTM_Converter(outputDir,refLat,refLon,refHeight,left,right,bottom,top):
    import os
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import xarray as xr
    from scipy.interpolate import RectBivariateSpline, griddata
    import utm
    from stl import mesh  # install with `pip install numpy-stl`
    from terrain import SRTM
    # # Terrain-Resolved Domain Setup
    # This notebook will generate a surface geometry file to be used by the microscale solver (e.g., SOWFA) to conform the solver domain to the terrain (e.g., with `moveDynamicMesh`).
    # 
    # Notebook written by Eliot Quon, modified by Regis Thedin\
    # {eliot.quon,regis.thedin}@nrel.gov

    # ## 1. User input
    # Output directory (absolute dir)
    outdir = outputDir
    # ### 1.1 Mesoscale parameters
    # Get lower resolution WRF terrain data
    getWRFdata = False
    # WRF reference solution to blend from low-res (WRF) at inflow boundaries to full SRTM resolution
    wrfout = '/home/rthedin/MMC/WFIP2Region/20161121_YSU/wrfout_d01_2016-11-22_00:00:00'
    # ### 1.2 Microscale parameters
    # 
    # First, set the resolution of the data the STL will be created from, and the output resolution.
    product = 'SRTM1' # SRTM1 | SRTM3 (30- and 90-m DEM)
    ds = 10. # output resolution


    # The following cell should be modified by the user, following the examples given. The cell contains information about the actual location to be saved.
    # 
    # The `refloc` variable is the location corresponding to (0,0) in SOWFA.
    # 
    # Set the fringe on each side of the domain. This fringe region will used to blend the high-resolution SRTM domain data with either i) low-reosolution WRF (mesoscale) digital elevation model (DEM), or ii) flat.
    # 
    # If `getWRFdata` above is `True`, then blending to mesoscale will occur; otherwise, the domain will be blended to flat. If blending to flat, the user can specify an extra fringe region of completely flat domain (`fringe_flat`). Additionally, if blending to flat, the terrain surface data can be shifted vertically such that the flat region is at $z=0$ by setting `shiftFlatToZero` to `True`.
    # 
    # With respect to the bounding box, it is nice to have the boundaries exactly where the mesh would go because of the blending. For instance, a 5x5 km domain needs to match all the levels of cells: 20, 40, 80, 160 m. Essentially, 5000/160 needs to be an integer number, the same way 5000/80 needs to be as well. However, we only really need to match the coarsest resolution because they are multiple. Finally, an extra fringe of width `ds` (the resolution set above) is added, half on each side. The desired bounding box should go into the `xmin`,`xmax`,`ymin`,`ymax` variables, ignoring the `ds` addition. This extra fringe is to ensure that the STL is slightly larger than the bounding box that will be set on the microscale solver (needed in OpenFOAM to avoid numerical issues).


    # For WFIP2 region
    # - blend to flat
    #refloc = (45.638004, -120.642973, 495) # biglow PS12 met mast
    refloc=(refLat,refLon,refHeight)
    xmin,xmax = -left-ds/2,right+ds/2
    ymin,ymax = -bottom-ds/2,top+ds/2
    fringe_flat=150
    shiftFlatToZero=True
    fringe_w = 3000
    fringe_s = 3000
    fringe_n = 3000
    fringe_e = 3000
    case = f'wfip_xm{abs(int(xmin))}to{int(xmax)}_ym{abs(int(ymin))}to{int(ymax)}_blendFlat3N3S3E3W_ff{fringe_flat}'

    # - blend to WRF
    # refloc = (45.638004, -120.642973, 495) # biglow PS12 met mast
    # xmin,xmax = -15000-ds/2, 15720+ds/2
    # ymin,ymax = -5000-ds/2, 15160+ds/2
    # fringe_flat=0
    # shiftFlatToZero=False
    # fringe_w = 3000
    # fringe_s = 3000
    # fringe_n = 3000
    # fringe_e = 3000
    # case = f'wfip_xm{abs(int(xmin))}to{int(xmax)}_ym{abs(int(ymin))}to{int(ymax)}_blendWRF3N3S3E3W'


    # ## 2. Create output surface

    # In[11]:


    x1 = np.arange(xmin, xmax+ds, ds)
    y1 = np.arange(ymin, ymax+ds, ds)
    xsurf,ysurf = np.meshgrid(x1, y1, indexing='ij')


    # In[12]:


    print('The output bounding box is')
    print('xmin: ',xsurf[0,0], '\nxmax: ',xsurf[-1,-1])
    print('ymin: ',ysurf[0,0], '\nymax: ',ysurf[-1,-1])


    # ## 3. Get the high-resolution terrain

    # In[13]:


    # Terrain region to clip from the digital elevation model (DEM)
    srtm_bounds = west, south, east, north = (refloc[1]-0.5, refloc[0]-0.4, refloc[1]+0.62, refloc[0]+0.42)

    # this will be downloaded:
    srtm_output=f'{outdir}/{case}.tif' # need absolute path for GDAL


    # In[14]:


    srtm = SRTM(srtm_bounds, fpath=srtm_output, product=product)


    # In[15]:


    #get_ipython().run_line_magic('time', 'srtm.download()')
    # CPU times: user 3.53 ms, sys: 12.7 ms, total: 16.2 ms
    # Wall time: 8.74 s


    # In[16]:
    srtm.download()
    x,y,z = srtm.to_terrain()
    #get_ipython().run_cell_magic('time', '', "# original SRTM terrain stored as 'z'\nx,y,z = srtm.to_terrain()\n")


    # In[17]:


    # get reference location to use as origin
    xref,yref,_,_ = utm.from_latlon(*refloc[:2], force_zone_number=srtm.zone_number)


    # In[18]:


    vmin,vmax = 1500,2500

    fig,ax = plt.subplots(figsize=(12,8))
    cm = ax.pcolormesh(x-xref, y-yref, z, cmap='terrain')#,vmin=vmin,vmax=vmax)
    cb = fig.colorbar(cm,ax=ax)
    cb.set_label('elevation [m]',fontsize='x-large')
    ax.tick_params(labelsize='large')
    ax.set_xlabel('easting [m]')
    ax.set_ylabel('northing [m]')
    ax.set_title(f'{product} DEM projection')
    ax.axis('scaled')

    # bounding box for microscale region
    les = Rectangle((xmin,ymin), xmax-xmin, ymax-ymin, edgecolor='r', lw=3, facecolor='0.5', alpha=0.5)
    ax.add_patch(les)
    #plt.show()

    # ### 3.1 Downscale to output grid

    # In[19]:


    interpfun = RectBivariateSpline(x[:,0]-xref, y[0,:]-yref, z)


    # In[20]:


    # resampled SRTM data stored in 'zsrtm'
    zsrtm = interpfun(x1,y1,grid=True)



    # In[21]:


    fig,ax = plt.subplots(figsize=(12,8))
    cm = ax.pcolormesh(xsurf, ysurf, zsrtm, cmap='terrain')#,vmin=vmin,vmax=vmax)
    cb = fig.colorbar(cm,ax=ax)
    cb.set_label('elevation [m]',fontsize='x-large')
    ax.tick_params(labelsize='large')
    ax.set_xlabel('easting [m]')
    ax.set_ylabel('northing [m]')
    ax.set_title(f'{product} terrain height')
    ax.axis('scaled')

    fig.savefig(f'{outdir}/elevation_srtm_{case}.png',dpi=150,bbox_inches='tight')


    # ## 4. Get the low-resolution terrain from the mesoscale
    # This part is only relevant if the user chose to blen the high-resolution SRTM terrain data with WRF

    # In[22]:


    # Open the dataset
    if getWRFdata is True:
        if os.path.isfile(wrfout) is True:
            wrf = xr.open_dataset(wrfout)
            wrf['HGT']
        else:
            print('WRF input does not exist')
            sys.exit(1)


    # In[23]:


    if getWRFdata is True:
        # wrf fields
        hgt = wrf['HGT'][0,:,:]
        xlat = wrf.coords['XLAT'][0,:,:]
        xlon = wrf.coords['XLONG'][0,:,:]


    # In[24]:


    if getWRFdata is True:
        get_ipython().run_line_magic('time', 'output_lat, output_lon = srtm.to_latlon(xsurf.ravel()+xref, ysurf.ravel()+yref)')
        # CPU times: user 3.38 s, sys: 157 ms, total: 3.54 s
        # Wall time: 3.54 s


    # In[25]:


    if getWRFdata is True:
        # interpolate to wrf surface elevation based on lat/lon (stored as 'zwrf')
        xi = np.stack((output_lat.ravel(),output_lon.ravel()),axis=-1)
        points = np.stack((xlat.values.ravel(),xlon.values.ravel()),axis=-1)
        values = hgt.values.ravel()
        zi = griddata(points,values,xi)
        zwrf = zi.reshape(xsurf.shape)
        # CPU times: user 3.14 s, sys: 931 ms, total: 4.07 s
        # Wall time: 2.57 s


    # In[26]:


    if getWRFdata is True:
        fig,ax = plt.subplots(figsize=(12,8))
        cm = ax.pcolormesh(xsurf, ysurf, zwrf, cmap='terrain')#,vmin=vmin,vmax=vmax)
        cb = fig.colorbar(cm,ax=ax)
        cb.set_label('elevation [m]',fontsize='x-large')
        ax.tick_params(labelsize='large')
        ax.set_xlabel('easting [m]')
        ax.set_ylabel('northing [m]')
        ax.set_title('WRF terrain height')
        ax.axis('scaled')

        fig.savefig(f'elevation_wrf_{os.path.basename(wrfout)}.png',dpi=150,bbox_inches='tight')


    # ## 5. Blend surface definitions

    # In[27]:


    # check distance from west boundary
    blend_w = np.ones(xsurf.shape)
    if fringe_w > 0:
        blend_w = np.minimum(np.maximum((xsurf-xmin-fringe_flat)/fringe_w, 0), 1)


    # In[28]:


    # check distance from east boundary
    blend_e = np.ones(xsurf.shape)
    if fringe_e > 0:
        blend_e = np.minimum(np.maximum((xmax-xsurf-fringe_flat)/fringe_e, 0), 1)


    # In[29]:


    # check distance from south boundary
    blend_s = np.ones(xsurf.shape)
    if fringe_s > 0:
        blend_s = np.minimum(np.maximum((ysurf-ymin-fringe_flat)/fringe_s, 0), 1)


    # In[30]:


    # check distance from north boundary
    blend_n = np.ones(xsurf.shape)
    if fringe_n > 0:
        blend_n = np.minimum(np.maximum((ymax-ysurf-fringe_flat)/fringe_n, 0), 1)


    # In[31]:


    # combine blending functions
    blend = blend_w * blend_e * blend_s * blend_n


    # In[32]:


    fig,ax = plt.subplots(figsize=(12,8))
    cm = ax.pcolormesh(xsurf, ysurf, blend, cmap='magma')
    cb = fig.colorbar(cm,ax=ax)
    ax.tick_params(labelsize='large')
    ax.set_xlabel('easting [m]')
    ax.set_ylabel('northing [m]')
    ax.set_title('blending function')
    ax.axis('scaled')


    # In[33]:


    # create flat surface to be blended
    # SRTM data is unlikely to be around the z=0 mark, so get the average 
    z0 = np.amin(zsrtm) #0 #np.mean(zsrtm)
    zflat = np.full(zsrtm.shape,z0)


    # In[34]:


    # surface to blend
    if getWRFdata is True:
        zlowres = zwrf
    else:
        zlowres = zflat


    # In[35]:


    # now, blend the high/low resolution elevations
    zblend = blend*zsrtm + (1-blend)*zlowres


    # In[36]:


    fig,ax = plt.subplots(figsize=(12,8))
    cm = ax.pcolormesh(xsurf, ysurf, zblend, cmap='terrain')#,vmin=vmin,vmax=vmax)
    cb = fig.colorbar(cm,ax=ax)
    cb.set_label('elevation [m]',fontsize='x-large')
    ax.tick_params(labelsize='large')
    ax.set_xlabel('easting [m]')
    ax.set_ylabel('northing [m]')
    ax.set_title('blended terrain height')
    ax.axis('scaled')
    fig.savefig(f'{outdir}/elevation_blended_{case}.png',dpi=150,bbox_inches='tight')


    # # 6. Shift terrain
    # Shifts the terrain data so that the flat borders are at $z=0$

    # In[37]:


    if shiftFlatToZero:
        zTerrainRef=zblend[0,0]
        zblend = zblend - zblend[0,0]
        case = case + '_flatz0'


    # In[38]:


    if shiftFlatToZero:
        fig,ax = plt.subplots(figsize=(12,8))
        cm = ax.pcolormesh(xsurf, ysurf, zblend, cmap='terrain')#,vmin=vmin,vmax=vmax)
        cb = fig.colorbar(cm,ax=ax)
        cb.set_label('elevation [m]',fontsize='x-large')
        ax.tick_params(labelsize='large')
        ax.set_xlabel('easting [m]')
        ax.set_ylabel('northing [m]')
        ax.set_title('shifted terrain height')
        ax.axis('scaled')
        fig.savefig(f'{outdir}/elevation_blended_{case}.png',dpi=150,bbox_inches='tight')


    # ## 6. Write out terrain surface STL
    stlout = f'{outdir}/terrain.stl'
    # output 'zblend' surface - can skip blending step and just output 'zsrtm'
    Npts = np.prod(xsurf.shape)
    stlpoints = np.stack((xsurf.ravel(),
                        ysurf.ravel(),
                        zblend.ravel()),  # <-- output surface here
                        axis=-1)

    stlindices = np.reshape(np.arange(Npts), xsurf.shape)
    Nx,Ny = xsurf.shape
    Nfaces = (Nx-1)*(Ny-1)*2

    surf = mesh.Mesh(np.zeros(Nfaces, dtype=mesh.Mesh.dtype))
    iface = 0 
    for i in range(Nx-1):
        for j in range(Ny-1):
            surf.vectors[iface,0,:] = stlpoints[stlindices[i,j],:]
            surf.vectors[iface,1,:] = stlpoints[stlindices[i+1,j],:]
            surf.vectors[iface,2,:] = stlpoints[stlindices[i+1,j+1],:]
            surf.vectors[iface+1,0,:] = stlpoints[stlindices[i+1,j+1],:]
            surf.vectors[iface+1,1,:] = stlpoints[stlindices[i,j+1],:]
            surf.vectors[iface+1,2,:] = stlpoints[stlindices[i,j],:]
            iface += 2
    assert (iface == Nfaces)
    #get_ipython().run_cell_magic('time', '', 'Nx,Ny = xsurf.shape\nNfaces = (Nx-1)*(Ny-1)*2\n\nsurf = mesh.Mesh(np.zeros(Nfaces, dtype=mesh.Mesh.dtype))\n\n#\n# manually define triangular faces for this simple quad mesh\n#\n# for iface, f in enumerate(faces):\n#     for dim in range(3):\n#         surf.vectors[iface][dim] = vertices[f[dim],:]\niface = 0 \nfor i in range(Nx-1):\n    for j in range(Ny-1):\n        surf.vectors[iface,0,:] = stlpoints[stlindices[i,j],:]\n        surf.vectors[iface,1,:] = stlpoints[stlindices[i+1,j],:]\n        surf.vectors[iface,2,:] = stlpoints[stlindices[i+1,j+1],:]\n        surf.vectors[iface+1,0,:] = stlpoints[stlindices[i+1,j+1],:]\n        surf.vectors[iface+1,1,:] = stlpoints[stlindices[i,j+1],:]\n        surf.vectors[iface+1,2,:] = stlpoints[stlindices[i,j],:]\n        iface += 2\nassert (iface == Nfaces)\n# CPU times: user 27.5 s, sys: 182 ms, total: 27.7 s\n# Wall time: 27.7 s\n')
    dpath = os.path.dirname(stlout)
    if (not dpath == '') and (not os.path.isdir(dpath)):
        os.makedirs(dpath)
        print('Created',dpath)
        
    surf.save(stlout)
    # surf.save(stlout, mode=mesh.stl.ASCII) # if ASCII STL is needed
    print('Saved',stlout)
    print(zblend[0,0])
    return xref,yref,zTerrainRef,srtm






