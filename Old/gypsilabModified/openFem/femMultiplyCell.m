function M = femMultiplyCell(varargin)
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
%|    #    |   FILE       : femMultiplyCell.m                             |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Multiply finite elements cells                |
%|  `---'  |                                                              |
%+========================================================================+

if (nargin == 1)
    M = varargin{1};
    
elseif (nargin == 2)
    A = varargin{1};
    B = varargin{2};
    
    if iscell(A) && iscell(B)
        M = A{1} * B{1} + ...
            A{2} * B{2} + ...
            A{3} * B{3} ;
        
    elseif iscell(A) && ~iscell(B)
        M{1} = A{1} * B ;
        M{2} = A{2} * B ;
        M{3} = A{3} * B ;
        
    elseif ~iscell(A) && iscell(B)
        M{1} = A * B{1} ;
        M{2} = A * B{2} ;
        M{3} = A * B{3} ;
        
    elseif ~iscell(A) && ~iscell(B)
        M = A * B;
        
    else
        error('femMultiplyCell.m : unavailable case')
    end
    
elseif (nargin == 3)
    A = varargin{1};
    B = varargin{2};
    C = varargin{3};
    
    if iscell(A) && iscell(B) && iscell(C)
        M = A{1} * B{2} * C{3} - ...
            A{1} * B{3} * C{2} + ...
            A{2} * B{3} * C{1} - ...
            A{2} * B{1} * C{3} + ...
            A{3} * B{1} * C{2} - ...
            A{3} * B{2} * C{1} ;
        
    elseif iscell(A) && iscell(B) && ~iscell(C)
        M = ( A{1} * B{1} + ...
            A{2} * B{2} + ...
            A{3} * B{3} ) * C;
        
    elseif iscell(A) && ~iscell(B) && iscell(C)
        M = A{1} * B * C{1} + ...
            A{2} * B * C{2} + ...
            A{3} * B * C{3} ;
        
    elseif ~iscell(A) && iscell(B) && iscell(C)
        M = A * ( B{1} * C{1} + ...
            B{2} * C{2} + ...
            B{3} * C{3} );
        
    elseif ~iscell(A) && ~iscell(B) && iscell(C)
        M{1} = A * B * C{1};
        M{2} = A * B * C{2};
        M{3} = A * B * C{3};
        
    elseif ~iscell(A) && iscell(B) && ~iscell(C)
        M{1} = A * B{1} * C;
        M{2} = A * B{2} * C;
        M{3} = A * B{3} * C;
        
    elseif iscell(A) && ~iscell(B) && ~iscell(C)
        M{1} = A{1} * B * C;
        M{2} = A{2} * B * C;
        M{3} = A{3} * B * C;
        
    elseif ~iscell(A) && ~iscell(B) && ~iscell(C)
        M = A * B * C;
    else
        error('femMultiplyCell.m : unavailable case')
    end
    
else
    error('femMultiplyCell.m : unavailable case')
end
end
