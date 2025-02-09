function [G,M] = stokesletEQ(Gamma,phi,viscosity,alpha,tol)
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
%|    #    |   FILE       : stokesletEQ.m                                 |
%|    #    |   VERSION    : 0.10                                          |
%|   _#_   |   AUTHOR(S)  : Luca Berti                                    |
%|  ( # )  |   CREATION   : 25.12.2018                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+
% COMPUTE THE SINGLE-LAYER OPERATOR ASSOCIATED WITH THE STOKESLET
% Input -> Gamma: integration domain
%       -> phi: finite elements
%       -> viscosity: fluid viscosity
%       -> alpha: multiplicative factor for the normal \otimes normal
%          correction on the single layer operator
% Output -> G: single-layer matrix
%        -> M: mass matrix
if isempty(tol) == 1
    % Solution without H-matrices - no parallel computations
%     G = cell(3,3);
%     Ntensor = cell(3,3);
%     M = cell(3,3);
%     Mass = integral(Gamma,phi,phi);
%     ONE = @(X,Y) ones(size(X,1),1);
%     for i=1:3
%         for j=i:3
%             name = ['[ij/r+rirj/r^3]',num2str(i),num2str(j)];
%             green = @(X,Y) 1/(8*pi*viscosity)*femGreenKernel(X,Y,name,[]);
%             G{i,j} = integral(Gamma,Gamma,phi,green,phi);
%             G{i,j} = G{i,j} + 1/(8*pi*viscosity)*regularize(Gamma,Gamma,phi,name,phi);            
%             Ntensor{i,j} = integral(Gamma,Gamma,ntimes(phi,i),ONE,ntimes(phi,j));
%             M{i,j} = sparse(size(Mass,1),size(Mass,2));
%             G{j,i} = G{i,j};
%             M{j,i} = M{i,j};
%             Ntensor{j,i} = Ntensor{i,j};
%         end
%         M{i,i} = Mass;
%     end
%     G = cell2mat(G);
%     Ntensor=cell2mat(Ntensor);
%     G = G + alpha*Ntensor;
%     M = cell2mat(M);
    % Solution without H-matrices - parallel computations
    G = cell(6,1);
    M = cell(6,1);
    Ntensor = cell(6,1);
    ind    = { [1 1] ; [2 2] ; [3 3] ; [1 2] ; [1 3] ; [2 3] };
    Mass = integral(Gamma,phi,phi);
    ONE = @(X,Y) ones(size(X,1),1);
    for n = 1:6
        % Indices
        i = ind{n}(1);
        j = ind{n}(2);
        
        % Single layer : G = \int_gamma \int_gamma  1/(8pi) (\delta_ij/r + r_i*r_j/|r|^3)
        name  = ['[ij/r+rirj/r^3]',num2str(i),num2str(j)];
        green = @(X,Y) 1/(8*pi*viscosity) .* femGreenKernel(X,Y,name,[]);
        G{n}  = integral(Gamma,Gamma,phi,green,phi);
        Ntensor{n} = integral(Gamma,Gamma,ntimes(phi,i),ONE,ntimes(phi,j));
        % Regularization
        G{n} = G{n} + 1/(8*pi*viscosity) .* regularize(Gamma,Gamma,phi,name,phi);
        if (i == j)
            M{n} = Mass;
        else
            M{n} = sparse(size(Mass,1),size(Mass,2));
        end
    end    
    G = [G{1} , G{4} , G{5} ; G{4} , G{2} , G{6} ; G{5} , G{6} , G{3} ];
    M = [M{1} , M{4} , M{5} ; M{4} , M{2} , M{6} ; M{5} , M{6} , M{3} ];
    Ntensor = [Ntensor{1} , Ntensor{4} , Ntensor{5} ; Ntensor{4} , Ntensor{2} , Ntensor{6} ; Ntensor{5} , Ntensor{6} , Ntensor{3} ];
    %G = cell2mat(G);
    %Ntensor=cell2mat(Ntensor);
    G = G + alpha*Ntensor;
    %M = cell2mat(M);
else
    % Solution with H-matrices
    G = cell(6,1);
    M = cell(6,1);
    ind    = { [1 1] ; [2 2] ; [3 3] ; [1 2] ; [1 3] ; [2 3] };
    Mass = integral(Gamma,phi,phi);
    for n = 1:6
        % Indices
        i = ind{n}(1);
        j = ind{n}(2);
        
        % Single layer : G = \int_gamma \int_gamma  1/(8pi) (\delta_ij/r + r_i*r_j/|r|^3)
        name  = ['[ij/r+rirj/r^3]',num2str(i),num2str(j)];
        green = @(X,Y) 1/(8*pi*viscosity) .* femGreenKernel(X,Y,name,[]);
        G{n}  = integral(Gamma,Gamma,phi,green,phi,tol);
        
        % Regularization
        G{n} = G{n} + 1/(8*pi*viscosity) .* regularize(Gamma,Gamma,phi,name,phi);
        if (i == j)
            M{n} = Mass;
        else
            M{n} = sparse(size(Mass,1),size(Mass,2));
        end
    end
    G = [G{1} , G{4} , G{5} ; G{4} , G{2} , G{6} ; G{5} , G{6} , G{3} ];
    M = [M{1} , M{4} , M{5} ; M{4} , M{2} , M{6} ; M{5} , M{6} , M{3} ];
end





end