function Mh = hmxFusion(Mh)
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
%|    #    |   FILE       : hmxFusion.m                                   |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Fusion and recompress for full, sparse and    |
%|  `---'  |                low-rank leaves                               |
%+========================================================================+

%%% H-Matrix
if (Mh.typ == 0)
    % Check and convert leaf data
    typ = zeros(1,4);
    for i = 1:4
        % Children data        
        typ(i) = Mh.chd{i}.typ;
        dat    = Mh.chd{i}.dat;        
        tol    = Mh.chd{i}.tol;
        dim    = size(Mh.chd{i});

        % Full or sparse matrix
        if (typ(i) == 2)
            % Non zeros terms
            n = nnz(dat);
            
            % Sparse
            if issparse(dat)
                % Compress
                if (n == 0)
                    A    = zeros(dim(1),0);
                    B    = zeros(0,dim(2));
                    flag = 1;
                    
                % Full
                elseif (n > 1/4*prod(dim))
                    dat  = full(dat);
                    flag = 2;
                    
                % Sparse
                else
                    flag = 0;
                end
                
            % Full
            else
                % Compress
                if (n == 0)
                    A    = zeros(dim(1),0);
                    B    = zeros(0,dim(2));
                    flag = 1;
                    
                % Sparse
                elseif (n <= 1/4*prod(dim))
                    dat  = sparse(double(dat));
                    flag = 2;
                 
                % Compress
                else
                    [A,B,flag] = hmxSVD(dat,tol);
                end
            end
            
            % Update
            if (flag == 1)
                Mh.chd{i}.dat = {A,B};
                Mh.chd{i}.typ = 1;
                typ(i)        = 1;
            
            elseif (flag == 2)
                Mh.chd{i}.dat = dat;
            end
        end
    end
    
    % Low-rank fusion (QRSVD)
    if (sum(typ==1) == 4)
        % Rank for each leaf
        rk = zeros(1,4);
        nk = 0;
        for i = 1:4
            rk(i) = size(Mh.chd{i}.dat{1},2);
            nk    = nk + sum(size(Mh.chd{i}))*rk(i); 
        end
        
        % Low rank matrix
        A = zeros(size(Mh,1),sum(rk),class(Mh.pos{1}));
        B = zeros(sum(rk),size(Mh,2),class(Mh.pos{2}));
        j = 0;
        for i = 1:4
            A(Mh.row{i},j+(1:rk(i))) = Mh.chd{i}.dat{1};
            B(j+(1:rk(i)),Mh.col{i}) = Mh.chd{i}.dat{2};            
            j = j + rk(i);
        end
        
        % Recompression
        [A,B] = hmxQRSVD(A,B,Mh.tol);
        
        % Update
        if (sum(size(Mh))*size(A,2) <= nk)  
            Mh.typ = 1;
            Mh.row = cell(1,4);
            Mh.col = cell(1,4); 
            Mh.chd = cell(1,4);
            Mh.dat = {A,B};
        end
    end
    
%%% Unavalaible case 
else
    error('hmxFusion.m : unavailable case')
end
end
