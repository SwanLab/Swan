function Gxy = femGreenKernel(X,Y,green,k)
%+========================================================================+
%|                                                                        |
%|              OPENFEM - LIBRARY FOR FINITE ELEMENT METHOD               |
%|           openFem is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2017-2018.          |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : femGreenKernel.m                              |
%|    #    |   VERSION    : 0.55                                          |
%|   _#_   |   AUTHOR(S)  : M. Aussal & F. Alouges & M. Averseng          |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.05.2019                                    |
%| ( === ) |   SYNOPSIS   : Usefull green kernel functions                |
%|  `---'  |                                                              |
%+========================================================================+

% Security
if (size(X,2) ~= 3) || (size(Y,2) ~= 3)
    error('femGreenKernel.m : unavailable case')
end
if isempty(k)
    k = 0;
end

% Distances between particles
Rxy = sqrt( ...
    (X(:,1) - Y(:,1)).^2 + ...
    (X(:,2) - Y(:,2)).^2 + ...
    (X(:,3) - Y(:,3)).^2 );

% For empty wave-number
if isempty(k)
    k = 0;
end

% Green kernel definition
if strcmp(green,'[1/r]')
    Gxy = 1./Rxy;   
    
elseif strcmp(green,'dx[1/r]')
    Gxy = - 1 ./ (Rxy.^3);
    
elseif strcmp(green,'dy[1/r]')
    Gxy = 1 ./ (Rxy.^3);
    
elseif strcmp(green,'[exp(ikr)/r]')
    Gxy = exp(1i*k*Rxy)./Rxy;          
    
elseif strcmp(green,'dx[exp(ikr)/r]')
    Gxy = (1i*k - 1./Rxy) .* exp(1i*k.*Rxy) ./ (Rxy.^2);

elseif strcmp(green,'dy[exp(ikr)/r]')
    Gxy = - (1i*k - 1./Rxy) .* exp(1i*k.*Rxy) ./ (Rxy.^2);

elseif strcmp(green(1:end-1),'gradx[1/r]')    
    j = str2double(green(end));
    Gxy = - (X(:,j)-Y(:,j)) ./ (Rxy.^3);
    
elseif strcmp(green(1:end-1),'grady[1/r]')     
    j = str2double(green(end));
    Gxy = (X(:,j)-Y(:,j)) ./ (Rxy.^3);    
    
elseif strcmp(green(1:end-1),'gradx[exp(ikr)/r]')  
    j = str2double(green(end));
    Gxy = (1i*k - 1./Rxy) .* exp(1i*k.*Rxy) .* ...
        (X(:,j)-Y(:,j)) ./ (Rxy.^2);
    
elseif strcmp(green(1:end-1),'grady[exp(ikr)/r]')    
    j = str2double(green(end));
    Gxy = - (1i*k - 1./Rxy) .* exp(1i*k.*Rxy) .* ...
        (X(:,j)-Y(:,j)) ./ (Rxy.^2);
    
elseif strcmp(green,'[log(r)]')
    Gxy = log(Rxy);
    
elseif strcmp(green,'[H0(kr)]')
    Gxy = besselh(0,k*Rxy);
    
elseif strcmp(green(1:end-1),'gradx[log(r)]')
    j = str2double(green(end));
    Gxy = (X(:,j)-Y(:,j)) ./ (Rxy.^2);    
    
elseif strcmp(green(1:end-1),'grady[log(r)]')
    j = str2double(green(end));
    Gxy = - (X(:,j)-Y(:,j)) ./ (Rxy.^2);    
    
elseif strcmp(green(1:end-1),'gradx[H0(kr)]')
    j   = str2double(green(end));
    Gxy = - k * besselh(1,k*Rxy) .* (X(:,j)-Y(:,j)) ./ Rxy;
    
elseif strcmp(green(1:end-1),'grady[H0(kr)]')
    j   = str2double(green(end));
    Gxy = k * besselh(1,k*Rxy) .* (X(:,j)-Y(:,j)) ./ Rxy;
    
elseif strcmp(green(1:end-2),'[ij/r+rirj/r^3]')        
    i = str2double(green(end-1));
    j = str2double(green(end));
    Gxy = (i==j)./Rxy + (X(:,i)-Y(:,i)).*(X(:,j)-Y(:,j))./(Rxy.^3);
    
elseif strcmp(green(1:end-3),'[rirjrk/r^5]')      
    i = str2double(green(end-2));
    j = str2double(green(end-1));
    k = str2double(green(end));    
    Gxy = (X(:,i)-Y(:,i)).*(X(:,j)-Y(:,j)).*(X(:,k)-Y(:,k))./(Rxy.^5);
    
else
    error('Error in femGreenKernel.m : unknown green kernel')
end

% Singularity
if strcmp(green,'[exp(ikr)/r]')
    Gxy(Rxy<1e-12) = 0 + 1i*k;
elseif strcmp(green,'[H0(kr)]')
    gamma         = 0.5772156649;
    Gxy(Rxy<1e-12) = 1 + 1i*(2/pi*(gamma+log(k/2)));    
else
    Gxy(Rxy<1e-12) = 0;
end
end
