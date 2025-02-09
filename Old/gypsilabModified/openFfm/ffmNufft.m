function MV = ffmNufft(X,Y,V,iflag,tol)
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
%|    #    |   FILE       : ffmNufft.m                                    |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Preparation for Non Uniform FFT               |
%|  `---'  |                                                              |
%+========================================================================+

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
    % Prevention de coplanarite en X
    X(end+1,:)  = X(end,:) + tol;

    % Prevention de coplanarite en Y
    Y(end+1,:) = Y(end,:) + tol;
    V(end+1)   = 0;

    % Transformee de Fourier non uniforme
    if 1
        MV = nufft3d3(Y(:,1),Y(:,2),Y(:,3),V,iflag,tol,X(:,1),X(:,2),X(:,3));
    else
        opt.nthreads = 4;
        opt.fftw     = 0;
        MV = finufft3d3(Y(:,1),Y(:,2),Y(:,3),V,iflag,tol,X(:,1),X(:,2),X(:,3),opt);
    end
    
    % Supression dernier point
    MV = MV(1:end-1); 
end
end
