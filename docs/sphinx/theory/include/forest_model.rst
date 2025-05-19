Forest Model
--------------
The forest model provides an option to include the drag from forested regions to be included in the momentum equation. The 
drag force is calculated as follows: 

.. math::

   F_i= - C_d L(x,y,z) U_i | U_i |


Here :math:`C_d` is the coefficient of drag for the forested region and :math:`L(x,y,z)` is the leaf area density (LAD) for the 
forested region. A three-dimensional model for the LAD is usually unavailable and is also cumbersome to use if there are thousands
of trees. Two different models are available as an alternative: 

.. math::
   L=\frac{LAI}{h}

.. math:: 
   L(z)=L_m \left(\frac{h - z_m}{h - z}\right)^n  exp\left[n \left(1 -\frac{h - z_m}{h - z}\right )\right]

Here :math:`LAI` is the leaf area index and is available from measurements, :math:`h` is the height of the tree, :math:`z_m` is the location 
of the maximum LAD, :math:`L_m` is the maximum value of LAD at :math:`z_m` and :math:`n` is a model constant with values  6 (below :math:`z_m`) and 0.5 
(above :math:`z_m`), respectively. :math:`L_m` is computed by integrating the following equation: 

.. math::
   LAI = \int_{0}^{h} L(z) dz 

The simplified model with uniform LAD is recommended for forested regions with no knowledge of the individual trees. LAI values can be used from 
climate model look-up tables for different regions around the world if no local remote sensing data is available. 

