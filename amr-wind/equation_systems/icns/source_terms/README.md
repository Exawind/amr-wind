# AMR-Wind ICNS ABL Equation Source terms:
# Corilois & Geostrophic Forcing

## Coriolis Forcing

By definition the cross product between the coriolis parameter $\textbf{f}$
and the velocity vector $\textbf{u}$ is the coriolis acceleration 
(force) and is added to the ICNS eqiation to account for the rotation of 
the earth, moving from the inertial frame of reference to a rotating 
frame of reference.

$$ - \textbf{f} \times \textbf{u} $$

## Geostrophic Forcing

By definition the cross product between the coriolis parameter and
the geostrophic wind vector $\mathbf{U_g}$ is the mean
pressure gradient.

$$ - \nabla P = -\rho \textbf{f} \times \mathbf{U_g} $$

Where the coriolis parameter is the earths angular velocity
$\Omega$ projected onto any latitude $\theta$.

$$ \textbf{f} \equiv [0, 2\Omega\cos \theta, 2\Omega\sin \theta] $$


After taking the cross product of the coriolis parameter with either, the velocity
vector or the geostrophic wind vector, there can be three resulting components.
To align with what is typically done in the atmospheric community only the 
vertical component of the coriolis parameter required. So, a Boolean parameter
S is introduced. [ABL Inputs](https://exawind.github.io/amr-wind/user/inputs_ABL)
By default S = False to remove the vertical contributions to the forcings. However,
if S = True all three components are included in the forcings.