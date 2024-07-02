'''
This utility converts a STL file of the terrain into terrain.amrwind format. 
Currently it uses the STL triangulation to generate the points. 
In the future, the code can be modified to use python interpolaton tools. 
'''

import sys 
import pyvista as pv 
import numpy as np 

fileName=sys.argv[1]
mesh=pv.read(fileName)
print(mesh.bounds)
xmin=0
xmax=0
ymin=0
ymax=0
zmin=0
zmax=0
meshSize=16
dx=(xmax-xmin)
dy=(ymax-ymin)
dz=(zmax-zmin)
print(dx,dy,dz)
print(dx/meshSize,dy/meshSize,dz/meshSize)
target=open("terrain.amrwind",'w')
print(np.shape(mesh.points))
for i in range(0,(np.shape(mesh.points))[0]):
    if((mesh.points[i][2]-zmin)<20):
        print(mesh.points[i][2]-zmin)
    target.write("%g %g %g\n"%(mesh.points[i][0]-xmin,mesh.points[i][1]-ymin,mesh.points[i][2]-zmin))
target.close()
data=np.genfromtxt("terrain.amrwind")
mesh=pv.PolyData(data)
mesh['elevation']=data[:,2]
surf = mesh.delaunay_2d()
surf.save("terrain.vtk")
# mesh=pv.PolyData("terrain.vtk")
# tempmesh=mesh
# bounds=tempmesh.bounds
# for i in range(0,1):
#     print("Run ",i)
#     bounds=tempmesh.bounds
#     xmin=bounds[0]
#     xmax=bounds[1]
#     ymin=bounds[2]
#     ymax=bounds[2]
#     smoothData=tempmesh.smooth(108)    
#     for k in range(0,len(mesh.points)):
#         xc=smoothData.points[i][0]
#         yc=smoothData.points[i][1]
#         zc=smoothData.points[i][2]
#         if(xc>0.25*xmin and xc<0.75*xmax and yc>0.25*ymin and yc<0.75*ymax):
#             smoothData.points[i]=mesh.points[i]
#     tempmesh=smoothData
# print("Bounds:",bounds)
# smoothData.save("smoothTerrain.vtk")
# #smoothData=mesh
# target=open("terrain.amrwind","w")
# for i in range(0,data.shape[0]):
#     target.write("%g %g %g\n"%(smoothData.points[i,0],smoothData.points[i,1],smoothData.points[i,2]-np.amin(smoothData.points[i,2])))
# target.close()
