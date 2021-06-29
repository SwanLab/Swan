function Gxy = ffmGreenKernel(X,Y,green,k)
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
%|    #    |   FILE       : ffmGreenKernel.m                              |
%|    #    |   VERSION    : 0.6                                           |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Green kernel computation using string         |
%|  `---'  |                definition                                    |
%+========================================================================+

% Security
if (size(X,2) ~= 3) || (size(Y,2) ~= 3)
    error('ffmGreenKernel.m : unavailable case')
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
    
elseif strcmp(green,'[exp(ikr)/r]')
    Gxy = exp(1i*k*Rxy)./Rxy;          

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
    error('Error in ffmGreenKernel.m : unknown green kernel')
end

% Singularity
if strcmp(green,'[exp(ikr)/r]')
    Gxy(Rxy<1e-8) = 0 + 1i*k;
else
    Gxy(Rxy<1e-8) = 0;
end
end
