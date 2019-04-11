## SWAN LOGO (WIP) ##


## SWAN ## 

*"It is not enough for code to work."
― Robert C. Martin*

_Swan_ is a Topology Optimization software developed in Matlab. Currently, it supports the resolution of 2D and 3D topology optimization problems with density or level-set as design variables.
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

Intermediate material properties are determined by the following implemented interpolation schemes:
- SIMP
- SIMP-ALL



Swan's modular design allows to combine several functions to define different optimization problems. Each function can be used as a cost or a constraint in the optimization problem.  The functions that are currently implemented are:
- Compliance
- Volume
- Perimeter
- Non-self adjoint compliance (used to minimize/maximize displacements)
- Homogenized elasticity matrix (used in micro scale problems)

In terms of optimization techniques the following optimizers are implemented:
 - Density: Projected Gradient, MMA and IPOPT. 
 - Level-set: SLERP, Projected SLERP and Hamilton-Jacobi


## Contact ##

For any inquiries, please contact us by opening an issue.

Current active developers are: Àlex Ferrer (@FerrerFerreAlex), Marc Núñez (@marcnunezc) and Oriol Trujillo (@Trujillo94)

Previous developers are acknowledged: Ferran De la Fuente, Nacho Izquierdo Pérez, Raül Rubio Serrano, Albert Torres Rubio

## References 
PENDING (if needed): reference to MMA, IPOPT, P. SLERP, HJ

Bendsøe, M.P. Structural Optimization (1989) 1: 193. https://doi.org/10.1007/BF0165094

S. Amstutz, C. Dapogny, and A. Ferrer, *SIMP-ALL: a generalized SIMP method based
on topological and shape derivatives.* WCSMO12, 2017.

S. Amstutz and H. Andra, *A new algorithm for topology optimization using a level-set
method,* 2066.

