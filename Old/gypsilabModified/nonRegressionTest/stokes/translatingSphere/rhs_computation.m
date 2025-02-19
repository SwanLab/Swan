function U = rhs_computation(M,sphere_translational_velocity)
%+========================================================================+
%|                                                                        |
%|           This function uses the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT :                                                            |
%| PROPERTY  :                                                            |
%| LICENCE   :                                                            |
%| CONTACT   :                                                            |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : rhs_computation.m                             |
%|    #    |   VERSION    : 0.10                                          |
%|   _#_   |   AUTHOR(S)  : Luca Berti                                    |
%|  ( # )  |   CREATION   : 25.12.2018                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+
% COMPUTE THE RIGHT-HAND SIDE OF THE STOKES SINGLE LAYER PROBLEM
% ONLY P1 CASE TESTED FOR THE MOMENT

N_unknown = size(M,1)/3;

U = M*[sphere_translational_velocity(1)*ones(N_unknown,1);
     sphere_translational_velocity(2)*ones(N_unknown,1);
     sphere_translational_velocity(3)*ones(N_unknown,1)];

end