import numpy as np

x=np.linspace(0,1000,201)
y=np.linspace(0,1000,201)

target=open("terrain.amrwind","w")
for i in range(0,len(x)):
    for j in range(0,len(y)):
        if(x[i]>=400 and x[i]<=600 and y[j]>=400 and y[j]<=600):
            z=200
        else:
            z=0
        target.write("%g %g %g\n"%(x[i],y[j],z))
target.close()
        
