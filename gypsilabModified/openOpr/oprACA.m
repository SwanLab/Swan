function [A,B,flag] = oprACA(opr,k,X,Xnrm,Y,Ynrm,tol)
%+========================================================================+
%|                                                                        |
%|            OPENOPR - LIBRARY FOR SPECIFIC OPERATORS IN BEM             |
%|           openOpr is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2019.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : oprACA.m                                      |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Adaptative Cross Approximation, partial       |
%|  `---'  |                & total pivoting                              |
%+========================================================================+

% Partial pivoting with handle function
row = @(i) oprGreenKernel(opr,k,X(i,:),Xnrm(i,:),Y,Ynrm);
col = @(j) oprGreenKernel(opr,k,X,Xnrm,Y(j,:),Ynrm(j,:));

% First row (row index of maximum reference value)
i = 1;
B = row(i);

% First pivot (maximum of the first row) 
[~,j] = max(abs(B));
delta = B(j);
if (abs(delta) < 1e-12)
    delta = 1e-12;
end

% First column
A = col(j)./delta;

% Matrix dimensions
Nx = size(A,1);
Ny = size(B,1);

% Unused row indice
Ir = (1:Nx)';
Ir = [Ir(1:i-1);Ir(i+1:end)];

% Unused column indice
Ic = (1:Ny)';
Ic = [Ic(1:j-1);Ic(j+1:end)];

% Frobenius norm for the initial tensor product
An2   = (A'*A);
Bn2   = (B'*B);
Rn2   = An2 * Bn2;

% Iterative construction using frobenius and infinite error on reference
n = 1;
while (sqrt(An2)*sqrt(Bn2) > tol*sqrt(Rn2))
    % New pivot by maximum 
    [~,i] = max(abs(A(Ir,n)));
    new   = row(Ir(i)) - B*A(Ir(i),1:n).';
    [~,j] = max(abs(new(Ic)));
    delta = new(Ic(j));
    
    % New pivot by minimum 
    [~,iMin] = min(abs(A(Ir,n)));
    newMin   = row(Ir(iMin)) - B*A(Ir(iMin),1:n).';
    [~,jMin] = max(abs(newMin(Ic)));
    deltaMin = new(Ic(jMin));
    
    % Best pivot
    if (abs(deltaMin)>abs(delta))
        i     = iMin;
        new   = newMin;
        j     = jMin;
        delta = deltaMin;
    end
        
    % Security for delta
    if (abs(delta) < 1e-12)
        delta = 1e-12;
    end
    
    % Update row
    B(:,n+1) = new;
    
    % New column
    A(:,n+1) = (col(Ic(j)) - A*B(Ic(j),1:n).') ./ delta;
    
    % Update row and column indices
    Ir = [Ir(1:i-1);Ir(i+1:end)];
    Ic = [Ic(1:j-1);Ic(j+1:end)];
    
    % Incrementation
    n = n + 1;

    % Relative Frobenius error by block
    An2 = A(:,n)'*A(:,n);
    Bn2 = B(:,n)'*B(:,n);
    u   = B(:,n)'*B(:,1:n-1);
    v   = A(:,n)'*A(:,1:n-1);
    AB  = v*u.';
    Rn2 = Rn2 + 2*real(AB) + An2*Bn2;
    
    % Compression failed
    if (n*(Nx+Ny) > Nx*Ny)
        A    = [];
        B    = [];
        flag = 0;
        return
    end    
%     sqrt(An2)*sqrt(Bn2)/sqrt(Rn2)
%     norm(A*B.' - A(:,1:n-1)*B(:,1:n-1).','fro')/norm(A*B.','fro')
%     norm(A(:,n)*B(:,n).','fro')/norm(A*B.','fro')    
%     norm(ref,'inf')/nrm
%     '================================'
end

% B transposition lead to  A * B
B    = B.';
flag = 1;
end
