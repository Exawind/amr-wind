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
    product = 'SRTM1' # SRTM1 | SRTM3 (30- and 90-m DEM)
    ds = 10. # output resolution

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
    x1 = np.arange(xmin, xmax+ds, ds)
    y1 = np.arange(ymin, ymax+ds, ds)
    xsurf,ysurf = np.meshgrid(x1, y1, indexing='ij')
    print('The output bounding box is')
    print('xmin: ',xsurf[0,0], '\nxmax: ',xsurf[-1,-1])
    print('ymin: ',ysurf[0,0], '\nymax: ',ysurf[-1,-1])
    srtm_bounds = west, south, east, north = (refloc[1]-0.5, refloc[0]-0.4, refloc[1]+0.62, refloc[0]+0.42)
    srtm_output=f'{outdir}/{case}.tif' # need absolute path for GDAL
    srtm = SRTM(srtm_bounds, fpath=srtm_output, product=product)
    srtm.download()
    x,y,z = srtm.to_terrain()
    xref,yref,_,_ = utm.from_latlon(*refloc[:2], force_zone_number=srtm.zone_number)
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
    interpfun = RectBivariateSpline(x[:,0]-xref, y[0,:]-yref, z)
    zsrtm = interpfun(x1,y1,grid=True)
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
    # check distance from west boundary
    blend_w = np.ones(xsurf.shape)
    if fringe_w > 0:
        blend_w = np.minimum(np.maximum((xsurf-xmin-fringe_flat)/fringe_w, 0), 1)
    blend_e = np.ones(xsurf.shape)
    if fringe_e > 0:
        blend_e = np.minimum(np.maximum((xmax-xsurf-fringe_flat)/fringe_e, 0), 1)
    blend_s = np.ones(xsurf.shape)
    if fringe_s > 0:
        blend_s = np.minimum(np.maximum((ysurf-ymin-fringe_flat)/fringe_s, 0), 1)
    blend_n = np.ones(xsurf.shape)
    if fringe_n > 0:
        blend_n = np.minimum(np.maximum((ymax-ysurf-fringe_flat)/fringe_n, 0), 1)
    blend = blend_w * blend_e * blend_s * blend_n
    fig,ax = plt.subplots(figsize=(12,8))
    cm = ax.pcolormesh(xsurf, ysurf, blend, cmap='magma')
    cb = fig.colorbar(cm,ax=ax)
    ax.tick_params(labelsize='large')
    ax.set_xlabel('easting [m]')
    ax.set_ylabel('northing [m]')
    ax.set_title('blending function')
    ax.axis('scaled')
    z0 = np.amin(zsrtm) #0 #np.mean(zsrtm)
    zflat = np.full(zsrtm.shape,z0)
    zlowres = zflat
    zblend = blend*zsrtm + (1-blend)*zlowres
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
    if shiftFlatToZero:
        zTerrainRef=zblend[0,0]
        zblend = zblend - zblend[0,0]
        case = case + '_flatz0'
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
    if (not dpath == '') and (not os.path.isdir(dpath)):
        os.makedirs(dpath)
        print('Created',dpath)
    surf.save(stlout)
    # surf.save(stlout, mode=mesh.stl.ASCII) # if ASCII STL is needed
    print('Saved',stlout)
    print(zblend[0,0])
    return xref,yref,zTerrainRef,srtm






