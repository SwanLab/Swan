function M = hmxSpy(varargin)
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
%|    #    |   FILE       : hmxSpy.m                                      |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Spy H-Matrix architecture                     |
%|  `---'  |                                                              |
%+========================================================================+

% Input data
Mh = varargin{1};

% Initialize visu block
if (Mh.typ > 0)
    M        = sparse(size(Mh,1),size(Mh,2));
    M(:,1)   = 1;
    M(:,end) = 1;
    M(1,:)   = 1;
    M(end,:) = 1;
end

%%%% H-Matrix (recursion)
if (Mh.typ == 0)
    A = hmxSpy(Mh.chd{1},[]);
    B = hmxSpy(Mh.chd{2},[]);
    C = hmxSpy(Mh.chd{3},[]);
    D = hmxSpy(Mh.chd{4},[]);
    M = [A,B;C,D];

%%% Compressed leaf
elseif (Mh.typ == 1)      
    rk = size(Mh.dat{1},2);
    if (rk > 0)
        M(:,rk) = 1;
        M(rk,:) = 1;
    else
        M = 3 .* M;
    end
       
%%% Full leaf
elseif (Mh.typ == 2)
    if issparse(Mh.dat)
        M = 4 .* (M + (abs(Mh.dat)>0));
    else
        M = 2 .* M;
    end
    
%%% Unknown type    
else
    error('hmxSpy.m : unavailable case')
end

% Graphical representation
if (length(varargin) == 1)
    spy(M==1,'b')
    hold on
    spy(M==2,'r')
    spy(M==3,'g')
    spy(M==4,'m')
    hold off
    title('openHmx : H-Matrix structure');
end
end
