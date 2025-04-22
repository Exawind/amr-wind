.. _terrainmodel:

Terrain Model
--------------

An immersed boundary forcing method (IBFM) is used to represent the terrain. In this method,
the effect of the terrain is modeled using a forcing term in the momentum and energy equation. 

The forcing term in the momentum equation is given by: 

.. math::

   F_i = - \beta C_d u_i | u_i | 

Here :math:`\beta` is the volume fraction of the cell covered by terrain, :math:`C_d` is a drag
term  and :math:`u_i` is the wind speed. Currently, the volume fraction is 
computed as a 0 or 1 using a simple nearest cell algorithm at each grid level. Future, updates 
will incorporate the partial terrain overlap using the EB capability in AMReX. The calculation 
of the drag coefficient term and the forcing term for the energy equation can be found in 
`Muñoz‐Esparza, Domingo, et al.  (JAMS 2020) <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002141>`_.

The original formulation is designed for low Reynolds number cases and does not include a 
method for applying a wall function. We propose the use of a forcing function to include 
the wall effects. 

First, compute the friction velocity from location k+1: 

.. math::

   u_*= |u_i[k+1]| \frac {\kappa}{\log [(z_{k+1}-z_k)/z0]}

The expected wind speed at cell k is computed as follows: 

.. math::

   |u_n|= \frac{u_*}{\kappa} \log [0.5 (z_{k+1}-z_k)/z0]

The methodology can be extended to include stability functions in a straight forward manner. The forcing 
term is computed as 

.. math::

   F_i= - \frac {|u[k]| \hat{c} - |u_n|\hat{l}} {\tau}



Here :math:`\hat{c}=(1,1,1)` is the existing normal vector from the grid and :math:`\hat{l}=(ux,uy,0)/|u_n|` is the value 
from the log law. The calculation of :math:`\hat{l}` and :math:`u_*` can be modified in the future align with the normal (following 
the orange arrow below).

.. image:: ./images/terrain_normal.png
   :align: center
   :width: 30%
