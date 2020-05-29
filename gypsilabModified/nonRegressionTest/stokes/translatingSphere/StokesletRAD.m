function G_eval = StokesletRAD(Gamma,phi,viscosity,mesh_eval)
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
%|    #    |   FILE       : Stokeslet_RAD.m                               |
%|    #    |   VERSION    : 0.10                                          |
%|   _#_   |   AUTHOR(S)  : Luca Berti                                    |
%|  ( # )  |   CREATION   : 25.12.2018                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+
% COMPUTE THE STOKESLET-RELATED RADIATION MATRIX
% Input -> Gamma: integration domain
%       -> phi: finite elements
%       -> viscosity: fluid viscosity
%       -> mesh_eval: mesh of evaluation points
% Output -> G_eval: single-layer radiation matrix

    X = mesh_eval.vtx;
    G_eval = cell(3,3);
    
    for i=1:3
        for j=1:3
            name = ['[ij/r+rirj/r^3]',num2str(i),num2str(j)];
            green = @(X,Y) 1/(8*pi*viscosity)*femGreenKernel(X,Y,name,[]);
            G_eval{i,j} = integral(X,Gamma,green,phi);
        end
    end
    
    G_eval = cell2mat(G_eval);
end