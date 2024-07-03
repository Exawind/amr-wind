from netCDF4 import Dataset
import numpy as np

data=Dataset("post_processing/sampling00000.nc")

velocityX=data['volume1']['velocityx']
velocityY=data['volume1']['velocityy']
velocityZ=data['volume1']['velocityz']
blank=data['volume1']['terrainBlank']
for i in range(0,velocityX.shape[0]-1):
    print(np.mean(velocityX[i][:]),np.mean(velocityY[i][:]),np.mean(velocityZ[i][:]))


