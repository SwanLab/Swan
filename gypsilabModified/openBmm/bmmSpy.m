function M = bmmSpy(varargin)
%+========================================================================+
%|                                                                        |
%|               OPENBMM - LIBRARY FOR BLOCK MATRIX ALGEBRA               |
%|           openBmm is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : bmmSpy.m                                      |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+

% Input data
Mb = varargin{1};

% Initialize visu block
[m,n] = size(Mb.blk);
M     = cell(m,n);

% For each block
for i = 1:m
    for j = 1:n
        % Get block
        Mij = Mb.blk{i,j};
        
        % Spy sparse matrix
        if isnumeric(Mij) && issparse(Mij)
            M{i,j} = 4*(abs(Mij) > 0);
        
        % Spy full matrix
        elseif isnumeric(Mij)
            M{i,j} = 2*(abs(Mij) > 0);
         
        % Spy H-Matrix    
        elseif isa(Mij,'hmx')
            M{i,j} = hmxSpy(Mb.blk{i,j},[]);

        % Spy FFM-Matrix
        elseif isa(Mij,'ffm')
            M{i,j}                 = sparse(size(Mij,1),size(Mij,2));
            M{i,j}(:,floor(end/2)) = 5;
            
        % Spy Block-Matrix
        elseif isa(Mij,'bmm')
            M{i,j} = bmmSpy(Mb.blk{i,j},[]);

        % Unknown type
        else 
           error('bmmSpy.m : unavailable case')  
        end
            
        % Add block limits
        M{i,j}(1:2:end,1:2)       = 10;
        M{i,j}(1:2:end,end-1:end) = 10;
        M{i,j}(1:2,1:2:end)       = 10;
        M{i,j}(end-1:end,1:2:end) = 10;
    end
end

% Matrix format
M = cell2mat(M);

% Graphical representation
if (nargin == 1)
    spy(M==1,'b')   % compressed
    hold on
    spy(M==2,'r')   % full
    spy(M==3,'g')   % empty
    spy(M==4,'m')   % sparse
    spy(M==5,'c')   % matrix-vector
    spy(M==10,'k')  % block limit
    hold off
    title('libBmm : Block-Matrix structure');
end
end
