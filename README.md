## SWAN LOGO (WIP) ##


## SWAN ## 

*"It is not enough for code to work."
― Robert C. Martin*

_Swan_ is a Topology Optimization software developed in Matlab. Currently, it supports the resoultion of 2D and 3D topology opitmization problems with density or level-set as design variables.
Swan aims to offer an adaptable framework that allows to implement new optimization techniques or functionalities fastly, thanks to its modular design.  



<p align="center">
  <img src="https://github.com/SwanLab/Examples/blob/master/Videos/Video_ImpCantileverHexahedra_Case_1_1_1_32.gif" alt="Solution" style="width: 600px;"/>
  <sub>Minimization of compliance subject to volume for a 3D cantilever beam, using Level Set methods.</sub>
</p>



## Current features ##
Swan currently supports the resolution of **macro** and **micro**\* scale problems  defined in terms of the following design variables:
- Density
- Level Set

<sub>*Note: 3D micro scale problems are not available yet. </sub>

Swan's modular design allows to combine several functions to define different optimization problems. Each function can be used as a cost or a constriant in the optimization problem.  The functions that are currently implemented are:
- Compliance
- Volume
- Perimeter
- Non-self adjoint compliance (used to minimize/maximize displacements)
- Homogenized elasticity matrix (used in micro scale problems)


## Contact ##

For any inquiries, please contact us by opening an issue.

Current active developers are: Àlex Ferrer (@FerrerFerreAlex), Marc Núñez (@marcnunezc) and Oriol Trujillo (@Trujillo94)

Previous developers are acknowledged: Ferran De la Fuente, Nacho Izquierdo Pérez, Raül Rubio Serrano, Albert Torres Rubio
