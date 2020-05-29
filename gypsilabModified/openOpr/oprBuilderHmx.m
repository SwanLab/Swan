function Mh = oprBuilderHmx(opr,k,Xunk,Mx,X,Xnrm,Yunk,My,Y,Ynrm,tol)
%+========================================================================+
%|                                                                        |
%|            OPENOPR - LIBRARY FOR SPECIFIC OPERATORS IN BEM             |
%|           openOpr is part of the GYPSILAB toolbox for Matlab           |
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
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : oprBuilderHmx.m                               |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2018                                    |
%| ( === ) |   SYNOPSIS   : Finite element builder with low-rank          |
%|  `---'  |                approximation for handle function             |
%+========================================================================+

% Initialize H-Matrix
Mh = hmx(Xunk,Yunk,tol);

% Admissibility
[isfar,Xdim,Ydim] = hmxFar(Mh);

% Compression for far distances
if isfar
    % ACA compression for quadrature matrix
    if strcmp(opr,'nxK')
        [A,B,flag] = oprACAv(opr,k,X,Xnrm,Y,Ynrm,tol);
    else
        [A,B,flag] = oprACA(opr,k,X,Xnrm,Y,Ynrm,tol);
    end

    % Finite element integration
    if flag
        % Integration 
        if strcmp(opr,'H') || strcmp(opr,'T')
            for n = 1:length(Mx)
                Mx{n} = Mx{n} * A;
                My{n} = B * My{n};
            end
            A = cell2mat(Mx);
            B = cell2mat(My');
        elseif strcmp(opr,'nxK')    
            A = [ Mx{1}*A(:,:,2) , - Mx{1}*A(:,:,3) , ...
                Mx{2}*A(:,:,3) , - Mx{2}*A(:,:,1) , ...
                Mx{3}*A(:,:,1) , - Mx{3}*A(:,:,2) ] ;
            B = [ B(:,:,2)*My{3} ; B(:,:,3)*My{2} ; ...
                B(:,:,3)*My{1} ; B(:,:,1)*My{3} ; ...
                B(:,:,1)*My{2} ; B(:,:,2)*My{1} ] ;
        else
            A = Mx * A;
            B = B * My;
        end
        
        % Recompression
        [A,B] = hmxQRSVD(A,B,tol);
    end
    
else
    flag = 0;
end

    
%%% Compression
if flag
    % Type
    Mh.typ = 1;
    
    % Low-rank
    Mh.dat = {A,B};


%%%% Full or sparse for smallest box (stopping criterion)
elseif (min(size(Mh)) < 100)
    % Type
    Mh.typ = 2;
    
    % Quadrature matrix
    if strcmp(opr,'nxK')
        Gxy = zeros(size(X,1),size(Y,1),3);
        for j = 1:size(Gxy,2)
            Gxy(:,j,:) = oprGreenKernel(opr,k,X,Xnrm,Y(j,:),Ynrm(j,:));
        end
    else
        Gxy = zeros(size(X,1),size(Y,1));
        for j = 1:size(Gxy,2)
            Gxy(:,j) = oprGreenKernel(opr,k,X,Xnrm,Y(j,:),Ynrm(j,:));
        end
    end
    
    % Finite element integration 
    if strcmp(opr,'H') || strcmp(opr,'T')
        Mh.dat = zeros(size(Mh));
        for n = 1:length(Mx)
            Mh.dat = Mh.dat + Mx{n} * Gxy * My{n};
        end
    elseif strcmp(opr,'nxK')
        Mh.dat = Mx{1} * Gxy(:,:,2) * My{3} - Mx{1} * Gxy(:,:,3) * My{2} + ...
            Mx{2} * Gxy(:,:,3) * My{1} - Mx{2} * Gxy(:,:,1) * My{3} + ...
            Mx{3} * Gxy(:,:,1) * My{2} - Mx{3} * Gxy(:,:,2) * My{1} ;
    else
        Mh.dat = Mx * Gxy * My;
    end
              
    
%%% H-Matrix (recursion)
else
    % Type
    Mh.typ = 0;
    
    % Subdivision for X
    [I1,I2] = hmxSubdivide(Xunk,Xdim);
    Mh.row  = {I1 , I1 , I2 , I2};
    
    % Subdivision for Y
    [I1,I2] = hmxSubdivide(Yunk,Ydim);
    Mh.col  = {I1 , I2 , I1 , I2};

    % H-Matrix (recursion)
    for i = 1:4
        % Dof indices
        Ir = Mh.row{i};
        Ic = Mh.col{i};
        
        % Fem matrix subdivision and quadratures points indices
        [Mxchd,Ix] = oprSubdivideCell(Mx,Ir,'left');
        [Mychd,Iy] = oprSubdivideCell(My,Ic,'right');

        % Recursion
        Mh.chd{i} = oprBuilderHmx(opr,k,Xunk(Ir,:),Mxchd,X(Ix,:),Xnrm(Ix,:),...
            Yunk(Ic,:),Mychd,Y(Iy,:),Ynrm(Iy,:),tol);
    end
    
    % Fusion
    Mh = hmxFusion(Mh);
end
end


%         % Add H-Vector
%         Mv(Ix) = Mv(Ix) + tmp
%
%
% % Corrective Matrix-vector product for Stokes Stresslet
% if strcmp(opr,('Tc'))
%     % H-Matrix (recursion)
%     if (Mh.typ == 0)
%         
%     % Compressed leaf
%     elseif (Mh.typ == 1)
%          TODO
% %         Mv = Mx\(A*sum(B,2));
%         
%     % Full leaf
%     elseif (Mh.typ == 2)
%         Mv = sum(Gxy*My,2);
%     end
% end
