function MV = ffmInteractionsSparse(X,Xbox,Y,Ybox,V,Ibox,green,k,edg,tol)
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
%|    #    |   FILE       : ffmInteractionsSparse.m                       |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Sparse product for low-frequency compressible |
%|  `---'  |                leaves                                        |
%+========================================================================+

% Initialisation du produit Matrice-Vecteur
MV = zeros(size(X,1),1,class(V));

% Quadrature des interpolations lagrangiennes
[Xq,ii,jj,kk,xq] = ffmQuadrature(green,k,edg,tol);

% Unicite des vecteurs de translation
XY        = Ybox.ctr(Ibox(:,2),:) - Xbox.ctr(Ibox(:,1),:);
[~,Il,It] = unique(floor(XY*1e6),'rows','stable');
Nt        = length(Il);

% Fonctions de transfert
Tx = cell(Nt,1);
Ty = cell(Nt,1);
for i = 1:Nt    
    [Tx{i},Ty{i}] = ffmTransfert(Xq,XY(Il(i),:),green,k,edg,tol);
end

% Interpolations en Y
Vy = cell(size(Ybox.ind));
for i = unique(Ibox(:,2)')
    % Boite centree en Y dans le cube unitaire d'interpolation
    iy = Ybox.ind{i};
    ny = length(iy);
    Ym = (2/edg) .* (Y(iy,:) - ones(ny,1)*Ybox.ctr(i,:));
    
    % Interpolation de Lagrange (Ym->Xq)
    Vy{i} = ffmInterp(ii,jj,kk,xq,Ym,V(iy),-1);    
end

% Translations
TVy = cell(size(Xbox.ind));
for i = 1:length(TVy)
    TVy{i} = 0;
end
for i = 1:size(Ibox,1)
    TVy{Ibox(i,1)} = TVy{Ibox(i,1)} + Tx{It(i)} * (Ty{It(i)} * Vy{Ibox(i,2)});
end

% Interpolations en X
for i = unique(Ibox(:,1)')
    % Boite centree en X dans le cube unitaire d'interpolation
    ix = Xbox.ind{i};
    nx = length(ix);
    Xm = (2/edg) .* (X(ix,:) - ones(nx,1)*Xbox.ctr(i,:));

    % Interpolation de Lagrange (Xq->Xm)
    MV(ix) = ffmInterp(ii,jj,kk,xq,Xm,TVy{i},+1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xq,ii,jj,kk,xq] = ffmQuadrature(green,k,edg,tol)
% Nuages de reference emetteurs X
x       = (-0.5*edg:edg/(10-1):0.5*edg)';
[x,y,z] = ndgrid(x,x,x);
Xref    = [x(:) y(:) z(:)];
Nref    = size(Xref,1);

% Nuage de reference transmetteur avec translation minimale
r0   = [2*edg 0 0];
Yref = [Xref(:,1)+r0(1) , Xref(:,2:3)];

% Xref dans le cube unitaire d'interpolation
a    = [-0.5 -0.5 -0.5]*edg;
b    = [0.5 0.5 0.5]*edg;
Xuni = ones(Nref,1)*(2./(b-a)) .* (Xref-ones(Nref,1)*(b+a)./2);

% Solution de reference
V     = ones(Nref,1) + 1i;
[I,J] = ndgrid(1:Nref,1:1:Nref);
ref   = reshape(ffmGreenKernel(Xref(I,:),Yref(J,:),green,k),Nref,Nref) * V;

% Initialisation
sol = 1e6;
nq  = 1;

% Boucle sur l'ordre de quadrature 
while norm(ref-sol)/norm(ref) > tol
    % Incrementation
    nq = nq + 1;
    
    % Indices de construction de l'interpolateur de Lagrange
    [ii,jj,kk] = ndgrid(1:nq,1:nq,1:nq);
    
    % Points de quadratures Tchebitchev (1D et 3D)
    xq = cos( (2*(nq:-1:1)-1)*pi/(2*nq))';
    Xq = [xq(ii(:)) xq(jj(:)) xq(kk(:))];
    
    % Vecteur de la translation
    [TA,TB] = ffmTransfert(Xq,r0,green,k,edg,tol);
    
    % Interpolation de Lagrange (Ym->Yq)
    Vy = ffmInterp(ii,jj,kk,xq,Xuni,V,-1);
    
    % Transfert (Yq->Xq)
    TVy = TA * (TB * Vy);
    
    % Interpolation de Lagrange (Xq->Xm)
    sol = ffmInterp(ii,jj,kk,xq,Xuni,TVy,+1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TA,TB] = ffmTransfert(Xq,xy,green,k,edg,tol)
% Dimensions
Nq = size(Xq,1);

% Interpolants en X
a   = [0 0 0];
b   = [1 1 1]*edg;
Txq = (ones(Nq,1)*(b-a)/2).*Xq + ones(Nq,1)*(b+a)/2;

% Interpolants en Y
a   = xy + a;
b   = xy + b;
Tyq = (ones(Nq,1)*(b-a)/2).*Xq + ones(Nq,1)*(b+a)/2;

% Noyaux de green avec compression aux interpolants
[TA,TB,flag] = hmxACA(Txq,Tyq,green,k,tol/10);

% Calcul direct si necessaire
if ~flag
    [I,J] = ndgrid(1:Nq,1:1:Nq);
    TA    = reshape(ffmGreenKernel(Txq(I,:),Tyq(J,:),green,k),Nq,Nq);
    TB    = 1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MV = ffmInterp(ii,jj,kk,xq,X,V,iflag)
% Dimensions
Nx = size(X,1);
nq = length(xq);

% Interpolateur de Lagrange par dimension d'espace
Aq = cell(1,nq);
for l = 1:nq^2
    if isempty(Aq{ii(l)})
        Aq{ii(l)} = ones(Nx,3,class(V));
    end
    if ii(l) ~= jj(l)
        Aq{ii(l)} = Aq{ii(l)} .* (X-xq(jj(l)))./(xq(ii(l))-xq(jj(l)));
    end
end

% Interpolation (X->Xq) 
if (iflag == -1)
    MV = zeros(nq^3,1,class(V));
    for l = 1:nq^3
        MV(l) = (Aq{ii(l)}(:,1) .* Aq{jj(l)}(:,2) .* Aq{kk(l)}(:,3)).' * V;
    end

% Interpolation (Xq->X)     
elseif (iflag == +1)
    MV = zeros(Nx,1,class(V));
    for l = 1:nq^3
        MV = MV + Aq{ii(l)}(:,1) .* Aq{jj(l)}(:,2) .* Aq{kk(l)}(:,3) .* V(l);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
