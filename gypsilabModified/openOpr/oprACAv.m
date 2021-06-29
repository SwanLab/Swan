function [A,B,flag] = oprACAv(opr,k,X,Xnrm,Y,Ynrm,tol)
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
%|    #    |   FILE       : oprACAv.m                                     |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Adaptative Cross Approximation, partial       |
%|  `---'  |                & total pivoting                              |
%+========================================================================+

% Matrix dimensions
Nx = size(X,1);
Ny = size(Y,1);

% Initialize indices
Ir = (1:Nx)';
Ic = (1:Ny)';

% Partial pivoting with handle function
row = @(i) oprGreenKernel(opr,k,X(i,:),Xnrm(i,:),Y,Ynrm);
col = @(j) oprGreenKernel(opr,k,X,Xnrm,Y(j,:),Ynrm(j,:));

% First row (row index of maximum reference value)
i   = 1;
B   = row(i);
dim = size(B,3);

% First pivot (maximum of the first row) 
[~,j] = max(max(abs(B),[],3));
delta = B(j,:,:);

% Securitys
delta(abs(delta)<1e-12) = 1e-12;

% First column
A = col(j)./delta;

% Unused indices1
Ir = [Ir(1:i-1);Ir(i+1:end)];
Ic = [Ic(1:j-1);Ic(j+1:end)];

% Frobenius norm for the initial tensor product
Rn2 = zeros(dim,1);
err = 0;
for l = 1:dim
    An2    = A(:,1,l)'*A(:,1,l);
    Bn2    = B(:,1,l)'*B(:,1,l);
    Rn2(l) = An2*Bn2;
    err    = max(err,sqrt(An2)*sqrt(Bn2)/sqrt(Rn2(l)));
end

% Iterative construction using frobenius norm
n = 1;
while (err > tol)
    % New pivot by maximum 
    [~,i] = max(max(abs(A(Ir,n,:)),[],3));
    new   = row(Ir(i));
    for l = 1:dim
        new(:,:,l) = new(:,:,l) - B(:,:,l)*A(Ir(i),1:n,l).';
    end
    [~,j] = max(max(abs(new(Ic,:,:)),[],3));
    delta = new(Ic(j),:,:);

    % New pivot by minimum 
    [~,iMin] = min(min(abs(A(Ir,n,:)),[],3));
    newMin   = row(Ir(iMin));
    for l = 1:dim
        newMin(:,:,l) = newMin(:,:,l) - B(:,:,l)*A(Ir(iMin),:,l).';
    end
    [~,jMin] = max(max(abs(newMin(Ic,:,:)),[],3));
    deltaMin = newMin(Ic(jMin),:,:);
    
    % Best pivot
    if (norm(deltaMin(:)) > norm(delta(:)))
        i     = iMin;
        new   = newMin;
        j     = jMin;
        delta = deltaMin;
    end
    
    % Security
    delta(abs(delta)<1e-12) = 1e-12;
    
    % Update row
    B(:,n+1,:) = new;
    
    % New column
    new = col(Ic(j));
    for l = 1:dim
        new(:,:,l) = (new(:,:,l) - A(:,:,l)*B(Ic(j),1:n,l).')./delta(l);    %B(:,:,l)*A(Ir(i),:,l).';
    end
    A(:,n+1,:) = new;
    
    % Update row and column indices
    Ir = [Ir(1:i-1);Ir(i+1:end)];
    Ic = [Ic(1:j-1);Ic(j+1:end)];
    
    % Incrementation
    n = n + 1;

    % Relative Frobenius error by block
    err = 0;
    for l = 1:dim
        An2    = A(:,n,l)'*A(:,n,l);
        Bn2    = B(:,n,l)'*B(:,n,l);
        u      = B(:,n,l)'*B(:,1:n-1,l);
        v      = A(:,n,l)'*A(:,1:n-1,l);
        AB     = v*u.';
        Rn2(l) = Rn2(l) + 2*real(AB) + An2*Bn2;
        err    = max(err,sqrt(An2)*sqrt(Bn2)/sqrt(Rn2(l)));
    end
    
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
B    = permute(B,[2 1,3]);
flag = 1;
end
