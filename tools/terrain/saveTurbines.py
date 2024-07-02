import pyvista as pv 
print("Plotting Turbines")
x1=self.terrainX1.flatten(order='F')
x2=self.terrainX2.flatten(order='F')
x3=self.terrainX3.flatten(order='F')
data=np.column_stack([x1,x2,x3])
for i in range(0,len(self.turbineX1)):
        print("Plotting turbine",i+1)
        newDisk=pv.Disc(center=(self.turbineX1[i],self.turbineX2[i],self.turbineX3[i]+80),inner=8,outer=50,normal=(1.0, 0.0,0.0), r_res=1, c_res=24)
        #pl.add_mesh(newDisk)
        if(i==0):
            globalBox=newDisk
        else:
            localBox=newDisk
            tempbox=globalBox.merge([localBox])
            globalBox=tempbox
for i in range(0,len(self.turbineX3)):
    self.turbineX3[i]+=80
data=np.column_stack([self.turbineX1,self.turbineX2,self.turbineX3])
mesh1=pv.PolyData(data)
pl.add_mesh(mesh1,render_points_as_spheres=True,point_size=10)      
tempFile=Path(self.caseParent,self.caseName,"terrainTurbine","turbines.vtk")  
pv.save_meshio(tempFile.as_posix(),globalBox)
pl.view_xy()
pl.show_axes()
pl.show()