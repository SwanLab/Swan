function MV = ffmInteractionsSparseHF(X,Xbox,Y,Ybox,V,Ibox,green,k,edg,tol)
%+========================================================================+
%|                                                                        |
%|         OPENFFM - LIBRARY FOR FAST AND FREE MEMORY CONVOLUTION         |
%|           openFfm is part of the GYPSILAB toolbox for Matlab           |
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
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : ffmInteractionsSparseHF.m                     |
%|    #    |   VERSION    : 0.6                                           |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Sparse product for high-frequency compressible|
%|  `---'  |                leaves                                        |
%+========================================================================+

% Initialisation du produit Matrice-Vecteur
MV = zeros(size(X,1),1,class(V));

% Quadrature spherique de Geggenbauer
[Xq,Wq,l] = ffmQuadratureHF(k,edg,tol);

% Unicite des vecteurs de translation
XY        = Ybox.ctr(Ibox(:,2),:) - Xbox.ctr(Ibox(:,1),:);
[~,Il,It] = unique(floor(XY*1e6),'rows','stable');
Nt        = length(Il);

% Fonctions de transfert
Txy = cell(Nt,1);
for i = 1:Nt
    Txy{i} = ffmTransfertHF(Xq,Wq,XY(Il(i),:),k,l);
end

% Convolution en Y
Vy = cell(size(Ybox.ind));
for i = unique(Ibox(:,2)')
    % Boite centree en Y
    iy = Ybox.ind{i};
    ny = length(iy);
    Ym = Y(iy,:) - ones(ny,1)*Ybox.ctr(i,:);
    
    % Transformee de Fourier inverse (Ym->Xq)
    Vy{i} = ffmInterpHF(Xq,Ym,V(iy),-1,tol);
    
    % Derivation du noyau en Y
    if strcmp(green(1:end-1),'grady[exp(ikr)/r]') 
        j     = str2double(green(end));
        Vy{i} = -1i*Xq(:,j) .* Vy{i};
    end    
end

% Translations
TVy = cell(size(Xbox.ind));
for i = 1:length(TVy)
    TVy{i} = 0;
end
for i = 1:size(Ibox,1)
    TVy{Ibox(i,1)} = TVy{Ibox(i,1)} + Txy{It(i)}.*Vy{Ibox(i,2)};
end

% Convolutions en X
for i = unique(Ibox(:,1)')
    % Boite centree en X
    ix = Xbox.ind{i};
    nx = length(ix);
    Xm = X(ix,:) - ones(nx,1)*Xbox.ctr(i,:);
    
    % Derivation du noyau en X
    if strcmp(green(1:end-1),'gradx[exp(ikr)/r]') 
        j      = str2double(green(end));
        TVy{i} = 1i*Xq(:,j) .* TVy{i};
    end      
    
    % Transformee de Fourier (Xq->X)
    MV(ix) = ffmInterpHF(Xm,Xq,TVy{i},+1,tol);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xq,Wq,l] = ffmQuadratureHF(k,edg,tol)
% Ordre harmonique (E. Darve)
l = floor( abs(k)*sqrt(3)*edg - log(tol) );

% Quadrature de Gauss sur [-1,1]
u     = 1:l;
u     = u ./ sqrt(4*u.^2 - 1);
[V,x] = eig( diag(u,-1) + diag(u,+1) );
[x,I] = sort(diag(x));
w     = 2*V(1,I)'.^2;

% Quadrature de Gauss sur [0,1] par transformation lineaire
a = 0;   b = 1;
x = 0.5*(b-a)*x + 0.5*(a+b);
w = 0.5*(b-a)*w;

% Quadrature de Gauss en elevation (phi)
phi = acos(2*x(end:-1:1)-1) - pi/2;
Wp  = 2*w(end:-1:1)';

% Quadrature reguliere en azimut (theta)
Ntheta = 2*(l+1);
theta  = 2*pi/Ntheta*(0:Ntheta-1)';
Wt     = 2*pi/Ntheta*ones(Ntheta,1);

% Quadrature spherique par produit tensoriel
[theta,phi] = ndgrid(theta,phi);
theta       = theta(:);
phi         = phi(:);
Wq          = Wt * Wp;
Wq          = Wq(:);

% Coordonnees carthesienne et produit par le nombre d'onde
[x,y,z] = sph2cart(theta,phi,1);
Xq      = k*[x,y,z];

% Securite
if (l>3)
    % Harmoniques spheriques jusqu'a l'ordre 3
    n   = 1;
    Ylm = zeros(length(theta),16);
    for ll = 0:3
        Plm = sqrt((2*ll+1)/(4*pi)) * legendre(ll,cos(pi/2-phi),'sch').';
        for m = -ll:ll
            if m<0
                Ylm(:,n) = Plm(:,abs(m)+1) .* sin(abs(m).*theta);
            else
                Ylm(:,n) = Plm(:,abs(m)+1) .* cos(abs(m).*theta);
            end
            n = n+1;
        end
    end
    
    % Test de l'integration spherique des harmoniques (produit scalaire)
    [m,n] = size(Ylm);
    if norm( (Ylm' * spdiags(Wq,0,m,m) * Ylm) - eye(n) , 'inf') > 1e-5
        error('ffmQuadratureHF.m - error 1');
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = ffmTransfertHF(Xq,ws,xy,k,l)
% Fonctions de hankel spherique 1ere espece
besselhs = @(n,z) sqrt(pi./(2.*z)) .* besselh(n+0.5,1,z);

% Distance
kr = k*norm(xy);

% cosinus(angle incidence)
x = Xq*(-xy/kr)';

% Initialisation de la recursion sur (besselhs(kr) * Pl(cos(phi))
P0 = ones(size(x));
P1 = x;
T  = besselhs(0,kr).*P0 + 3i*besselhs(1,kr).*P1;

% Construction recursive des polynomes de Legendre
for ll = 2:l
    P2 = 1/ll .* ( (2*(ll-1)+1).*x.*P1 - (ll-1).*P0 );
    T  = T + (2*ll+1) * 1i^ll * besselhs(ll,kr) .* P2;
    P0 = P1;
    P1 = P2;
end

% Operateur de Translation
T = 1i*abs(k)/(4*pi) .* ws .* T;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MV = ffmInterpHF(X,Y,V,iflag,tol)
% Dimensions
Nx   = size(X,1);
Ny   = size(Y,1);
Nmax = max(Nx,Ny);
Nmin = min(Nx,Ny);

% DFT non uniforme
if (Nmin<100) || (Nmin<log(Nmax))
    MV = exp(1i*iflag*X*Y') * V;
    
% FFT non uniforme    
else 
    % Conversion de type
    if isa(X,'single') || isa(Y,'single') || isa(V,'single')
        sgl = 1;
        X = double(X); Y = double(Y); V = double(V);
    else
        sgl = 0;
    end
    
    % Initialisation produit Matrice-Vecteur
    MV = zeros(Nx,1) + 1i*zeros(Nx,1);
    
    % Prevention de coplanarite en X
    Nx         = Nx + 1;
    X(end+1,:) = X(end,:) + tol;
    MV(end+1)  = 0;

    % Prevention de coplanarite en X
    Ny         = Ny + 1;
    Y(end+1,:) = Y(end,:) + tol;
    V(end+1)   = 0;
    
    % Transformee de Fourier
    ier    = 0;
    mex_id = 'nufft3d3f90(i int[x], i double[], i double[], i double[], i dcomplex[], i int[x], i double[x], i int[x], i double[], i double[], i double[], io dcomplex[], io int[x])';
    MV     = nufft3d(mex_id,Ny,Y(:,1),Y(:,2),Y(:,3),V, iflag, tol, ...
        Nx,X(:,1),X(:,2),X(:,3),MV, ier, 1, 1, 1, 1, 1);
    
    % Supression dernier point
    MV = MV(1:end-1);
    
    % Conversion de type
    if sgl
        MV = single(MV);
    end    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
