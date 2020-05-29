function [A,B] = hmxLowrankSub(Mh,I,J)
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
%|    #    |   FILE       : hmxLowrankSub.m                               |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   :  Extract low-rank submatrix from H-Matrix     |
%|  `---'  |                                                              |
%+========================================================================+

%%% H-Matrix (recursion)
if (Mh.typ == 0)
    % Initialisation
    A = zeros(length(I),0,class(Mh.row{1}));
    B = zeros(0,length(J),class(Mh.row{1}));

    % Recursion
    for i = 1:4
        % Row indices
        [Ia,Ib] = ismember(I,Mh.row{i});
        Ib      = Ib(Ia);

        % Colums indices
        [Ja,Jb] = ismember(J,Mh.col{i});
        Jb      = Jb(Ja);
        
        % Computation for non empty indices
        if ~isempty(Ib) && ~isempty(Jb)
            % Children evaluation
            [Ai,Bi] = hmxLowrankSub(Mh.chd{i},Ib,Jb);
            
            % Low-rank addition for A
            if (size(Ai,2) > 0)
                tmp       = zeros(length(I),size(Ai,2),class(Mh.row{1}));
                tmp(Ia,:) = Ai;
                A         = [A , tmp];
                
                % Low-rank addition for B
                tmp       = zeros(size(Bi,1),length(J),class(Mh.row{1}));
                tmp(:,Ja) = Bi;
                B         = [B ; tmp];
            end
        end
    end

    % Recompression
    [A,B] = hmxQRSVD(A,B,Mh.tol);

%%% Compressed leaf
elseif (Mh.typ == 1)
    A = Mh.dat{1}(I,:);
    B = Mh.dat{2}(:,J);

%%% Full leaf
elseif (Mh.typ == 2)
    % Low-Rank conversion
    A = full(Mh.dat(I,J));
    B = eye(length(J));
    
    % Recompression
    [A,B] = hmxQRSVD(A,B,Mh.tol);

%%% Others    
else
    error('hmxLowrankSub.m : unavailable case')
end
end
