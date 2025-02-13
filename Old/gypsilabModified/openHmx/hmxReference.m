function [Ix,Iy,V,flag] = hmxReference(Nx,Ny,mat,rkMax)
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
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : hmxReference.m                                |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 25.11.2018                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Cross method for non zeros block search to    |
%|  `---'  |                compute reference matrix indices and values   |
%+========================================================================+

% Initialize output
Mxy = zeros(Nx,3);

% First block is full matrix
Ir = (1:Nx)';
Ic = (1:Ny)';

% Block search
n = 1;
while ~isempty(Ir) && ~isempty(Ic) && (n<=log(min(Nx,Ny))) && (n<rkMax)   % (n*(Nx+Ny) > Nx*Ny)
    % New row and non empty columns indices
    i   = Ir(1);
    row = mat(i*ones(size(Ic)),Ic);
    jdx = Ic(abs(row)>1e-12);
    
    % New column and non empty row indices
    if ~isempty(jdx)
        j   = jdx(1);
        col = mat(Ir,j*ones(size(Ir)));
        idx = Ir(abs(col)>1e-12);
    else
        idx = i;
    end
    
    % Update block research indices
    Ir = setdiff(Ir,idx);
    Ic = setdiff(Ic,jdx);
    
    % Compute reference indices and values
    if isempty(jdx)
        Mxy(idx,:) = [idx 1 0];
    else
        ind        = ceil(length(jdx)*rand(length(idx),1));
        Mxy(idx,:) = [idx jdx(ind) mat(idx,jdx(ind))];
    end
    
    % Incrementation
    n = n + 1;
end

% Remaining rows
if ~isempty(Ir) 
    if isempty(Ic)
        Ic = (1:Ny)';
    end
    ind       = ceil(length(Ic)*rand(length(Ir),1));
    Mxy(Ir,:) = [Ir Ic(ind) mat(Ir,Ic(ind))];
end

% Output format
Ix   = Mxy(:,1);
Iy   = Mxy(:,2);
V    = Mxy(:,3);
flag = 1;

% Security
if (norm(Ix-(1:Nx)','inf')>1e-12)
    error('hmxReference.m : unavailable case')
end
end