function MV = ffmProduct(X,Y,V,green,k,tol)
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
%|    #    |   FILE       : ffmProduct.m                                  |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Fast convolution product memory free using a  |
%|  `---'  |                recursive octree                              |
%+========================================================================+

% Dimensions utiles
Nx = size(X,1);
Ny = size(Y,1);

% Securites
if (size(X,2)~=3) || (size(Y,2)~=3) || (size(V,1)~=Ny) || (size(V,2)~=1)
    error('ffmProduct.m : uncumpatible data sizes, please check input')
end
if isempty(k)
    k = 0;
end

% Produit Matrice-Vecteur
if strcmp(green,'[exp(-ikxy)]')
    MV = ffmNufft(k.*X,Y,V,-1,tol);
    return
elseif strcmp(green(1:end-1),'gradx[exp(-ikxy)]')  
    j  = str2double(green(end));
    MV = (-1i*k*X(:,j)) .* ffmNufft(k.*X,Y,V,-1,tol);
    return
else
    MV = zeros(Nx,1,class(V));
end

% Definition des boites
Xmin = min(X);
Ymin = min(Y);
edg  = max( [max(X)-Xmin , max(Y)-Ymin] );
edg  = 1.01*edg;

% Structure de donnee en X
Xbox.ctr    = Xmin + 0.5*edg;
Xbox.ind{1} = (1:Nx)';
Xbox.nbr    = Nx;

% Structure de donnee en Y
Ybox.ctr    = Ymin + 0.5*edg;
Ybox.ind{1} = (1:Ny)';
Ybox.nbr    = Ny;

% Interaction entre les boites
Iprt = [1 1];

% Fonctions de hankel spherique 1ere espece
besselhs = @(n,z) sqrt(pi./(2.*z)) .* besselh(n+0.5,1,z);

% Clustering hierarchique par recurion
while (~isempty(Iprt)) && (mean(Xbox.nbr) > 100) && (mean(Ybox.nbr) > 100)    
    % Subdivision des boites (incrementation)
    [Xbox,Ix] = ffmSubdivide(X,Xbox,edg,0);
    [Ybox,Iy] = ffmSubdivide(Y,Ybox,edg,0);
    edg = 0.5*edg;    

    % Propagation des interactions aux enfants avec tri (super important)
    Ichd = [];
    for i = 1:8
        idx    = Ix(Iprt(:,1),i);
        boolx = (idx>0);
        for j = 1:8
            jdx  = Iy(Iprt(:,2),j);
            bool = logical(boolx .* (jdx>0));
            Ichd = [Ichd ; idx(bool) jdx(bool)];
        end
    end
    Ichd = sortrows(Ichd);
    
    % Booleen des interactions lointaines
    Ifar = sqrt(sum( (Xbox.ctr(Ichd(:,1),:)-Ybox.ctr(Ichd(:,2),:)).^2 , 2));
    Ifar = (Ifar>=2*edg-1e-6);
    
    % Indices des interactions a compresser ou propager
    Iprt = Ichd(~Ifar,:);
    Ifar = Ichd(Ifar,:);
    
    % Interactions avec compression
    if ~isempty(Ifar)
        % Ordre harmonique
        l = floor( abs(k)*sqrt(3)*edg - log(tol) );

        % Compression Geggenbauer ou Lagrange
        if (abs(besselhs(l,2*k*edg)) < 1/tol) && ...
                sum(strcmp(green(1:end-1),{'[exp(ikr)/r','gradx[exp(ikr)/r]','grady[exp(ikr)/r]'}))
            MV = MV + ffmInteractionsSparseHF(X,Xbox,Y,Ybox,V,Ifar,green,k,edg,tol);
        else
            MV = MV + ffmInteractionsSparse(X,Xbox,Y,Ybox,V,Ifar,green,k,edg,tol);
        end
%         MV = MV + ffmInteractionsFull(X,Xbox,Y,Ybox,V,Ifar,green,k);
    end
end

% Interactions sans compressions
if ~isempty(Iprt)
    MV = MV + ffmInteractionsFull(X,Xbox,Y,Ybox,V,Iprt,green,k);
end
end
