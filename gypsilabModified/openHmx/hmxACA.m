function [A,B,flag] = hmxACA(varargin)
%+========================================================================+
%|                                                                        |
%|         OPENHMX - LIBRARY FOR H-MATRIX COMPRESSION AND ALGEBRA         |
%|           openHmx is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
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
%|    #    |   FILE       : hmxACA.m                                      |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Adaptative Cross Approximation, partial       |
%|  `---'  |                & total pivoting                              |
%+========================================================================+

% Total pivoting with full matrix
if (nargin == 2)
    M     = varargin{1};
    tol   = varargin{2};
    rkMax = 1e6;
    row   = @(i) full(M(i,:));
    col   = @(i) full(M(:,i));
    mat   = @(i,j) M(sub2ind(size(M),i,j));
    
% Total pivoting with full matrix and maximum rank
elseif (nargin == 3)
    M     = varargin{1};
    tol   = varargin{2};
    rkMax = varargin{3};
    row   = @(i) full(M(i,:));
    col   = @(i) full(M(:,i));
    mat   = @(i,j) M(sub2ind(size(M),i,j));
    
% Partial pivoting with handle function    
elseif (nargin == 4) 
    X     = varargin{1};
    Y     = varargin{2};
    green = varargin{3};
    tol   = varargin{4};
    rkMax = 1e6;
    row = @(i) green(X(i,:),Y).';
    col = @(i) green(X,Y(i,:));
    mat = @(i,j) green(X(i,:),Y(j,:));
    
% Partial pivoting with ffm function
elseif (nargin == 5) && ischar(varargin{3})
    X     = varargin{1};
    Y     = varargin{2};
    green = varargin{3};
    k     = varargin{4};
    tol   = varargin{5};
    rkMax = 1e6;
    row   = @(i) ffmGreenKernel(X(i,:),Y,green,k).';
    col   = @(i) ffmGreenKernel(X,Y(i,:),green,k);
    mat   = @(i,j) ffmGreenKernel(X(i,:),Y(j,:),green,k);
    
% Partial pivoting with handle function and maximum rank   
elseif (nargin == 5) 
    X     = varargin{1};
    Y     = varargin{2};
    green = varargin{3};
    tol   = varargin{4};
    rkMax = varargin{5};
    row = @(i) green(X(i,:),Y).';
    col = @(i) green(X,Y(i,:));
    mat = @(i,j) green(X(i,:),Y(j,:));
    
else
    error('hmxACA : unavailable case')
end

% Dimensions
if exist('M','var')
    [Nx,Ny] = size(M);
else
    Nx = size(X,1);
    Ny = size(Y,1);
end

% Initialize indices for row and columns pivot (dicreasing)
Ir = (1:Nx)';
Ic = (1:Ny)';

% Initialize reference indices and values (fixed)
[Ix,Iy,ref,flag] = hmxReference(Nx,Ny,mat,rkMax);
if flag
    nrm = norm(ref,2);
else
    A    = [];
    B    = [];
    flag = 0;
    return
end

% First row (row index of maximum reference value)
[~,i]  = max(abs(ref));
B      = row(i).';

% First pivot (maximum of the first row) 
[~,j] = max(abs(B));
delta = B(j);
if (abs(delta) < 1e-12)
    delta = 1e-12;
end

% First column
A = col(j)./delta;

% Update row indices
Ir = [Ir(1:i-1);Ir(i+1:end)];
Ic = [Ic(1:j-1);Ic(j+1:end)];

% Frobenius norm for the initial tensor product
An2   = (A'*A);
Bn2   = (B'*B);
Rn2   = An2 * Bn2;

% Update references values
ref = ref - A(Ix).*B(Iy);

% Iterative construction using frobenius and infinite error on reference
n = 1;
while (sqrt(An2)*sqrt(Bn2) > tol*sqrt(Rn2)) || (norm(ref,2) > tol*nrm)
    % Row index (maximum of column or reference)
    [a,i] = max(abs(A(Ir,n)));
    [b,k] = max(abs(ref(Ir)));
    if (b>a*delta)
        i = k;
    end
    
    % New row
    new = row(Ir(i)).' - B*A(Ir(i),1:n).';
    
    % New pivot
    [~,j] = max(abs(new(Ic)));
    delta = new(Ic(j));
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
    
    % Update reference
    ref = ref - A(Ix,n).*B(Iy,n);
    
    % Compression failed
    if (n*(Nx+Ny) > Nx*Ny) || (n>=rkMax) 
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
