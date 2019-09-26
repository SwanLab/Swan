## Swan - Topology Optimization Laboratory ##

<p align="center">
  <img src="https://github.com/SwanLab/Utilities/blob/master/swan_logo_test2.png" alt="drawing" width="200"/>
</p>

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

In terms of optimization techniques, the following optimizers are implemented:
 - Density: Projected Gradient, MMA and IPOPT. 
 - Level-set: SLERP, Projected SLERP and Hamilton-Jacobi
 
Unconstrained optimizers are combined with an Aumgmented Lagrangian to solve constrained problems. 


## Contact ##

For any inquiries, please contact us by opening an issue.

Current active developers are: Àlex Ferrer (@FerrerFerreAlex), Marc Núñez (@marcnunezc) and Oriol Trujillo (@Trujillo94)

Previous developers are acknowledged: Ferran De la Fuente, Nacho Izquierdo Pérez, Raül Rubio Serrano, Albert Torres Rubio

## References 

M.  P.  Bendsøe,  Optimal  shape  design  as  a  material  distributionproblem, Structural optimization,  vol.  1,  no.  4,  pp.  193–202,  Dec1989. [Online]. Available: https://doi.org/10.1007/BF01650949

S. Amstutz, H. Andrä, A new algorithm for topology optimization using a level-set method, J. Comput. Phys. 216 (2) (2006) 573–588.

K. Svanberg, The method of moving asymptotes – a new method for structural optimization, International Journal for Numerical Methods in Engineering, 1987, 24, 359-373.

S.  Amstutz,  C.  Dapogny,  and  A.  Ferrer,  A  consistent  relaxation of   optimal   design   problems   for   coupling   shape   and   topological derivatives, Numerische Mathematik, vol. 140, no. 1, pp. 35–94, Sep 2018. [Online]. Available: https://doi.org/10.1007/s00211-018-0964-4

G.  Allaire,  F.  Jouve,  A.-M.  Toader,  Structural  optimization  using  sensitivity  analysis  and  a  level-set  method,  J.  Comput.  Phys.  194  (1)(2004)  363–393,  http://dx.doi.org/10.1016/j.jcp.2003.09.032.

A. Wächter and L. T. Biegler, On the Implementation of a Primal-Dual Interior Point Filter Line Search Algorithm for Large-Scale Nonlinear Programming, Mathematical Programming 106(1), pp. 25-57, 2006
