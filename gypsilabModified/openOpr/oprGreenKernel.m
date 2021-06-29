function Gxy = oprGreenKernel(opr,k,X,Xnrm,Y,Ynrm)
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
%|    #    |   FILE       : oprGreenKernel.m                              |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Operator green kernel functions               |
%|  `---'  |                                                              |
%+========================================================================+

% Distances between particles
Rxy = sqrt( ...
    (X(:,1) - Y(:,1)).^2 + ...
    (X(:,2) - Y(:,2)).^2 + ...
    (X(:,3) - Y(:,3)).^2 );

%%% [exp(ikr)/r]
if strcmp(opr,'S') || strcmp(opr,'H') || strcmp(opr,'T')
    Gxy = exp(1i*k*Rxy)./Rxy; 
    Gxy(Rxy<1e-12) = 0 + 1i*k;
    

%%% dny[exp(ikr)/r]    
elseif strcmp(opr,'D')
    Gxy = - (1i*k - 1./Rxy) .* exp(1i*k.*Rxy) ./ (Rxy.^2) .* sum((X-Y).*Ynrm,2);
    Gxy(Rxy<1e-12) = 0;
    
    
%%% dnx[exp(ikr)/r]    
elseif strcmp(opr,'Dt')
    Gxy = (1i*k - 1./Rxy) .* exp(1i*k.*Rxy) ./ (Rxy.^2) .* sum((X-Y).*Xnrm,2);
    Gxy(Rxy<1e-12) = 0;
    
    
%%% grady[exp(ikr)/r]    
elseif strcmp(opr,'nxK')
    Gxy = - (1i*k - 1./Rxy) .* exp(1i*k.*Rxy) ./ (Rxy.^2);
    Gxy(Rxy<1e-12) = 0;
    Gxy = Gxy .* (X-Y);
    Gxy = permute(Gxy,[1 3 2]);

    
%%% [ij/r+rirj/r^3]
elseif strcmp(opr,'G')
    I = X(:,4);
    J = Y(:,4);
    if (length(I)==1) && (size(X,1)==1)
        jdx = [J==1;J==2;J==3];
        Gxy = (I==J)./Rxy + (X(I)-Y(:,I)) .* (X(J)'-Y(jdx)) ./ (Rxy.^3);
    elseif (length(J)==1) && (size(Y,1)==1)
        idx = [I==1;I==2;I==3];
        Gxy = (I==J)./Rxy + (X(idx)-Y(I)') .* (X(:,J)-Y(J)) ./ (Rxy.^3);
    else
        error('oprGreenKernel.m : unavailable case')
    end
    Gxy(Rxy<1e-12) = 0;
  
    
%%% [rirjr.n/r^5]
elseif strcmp(opr,'Ts')
    I = X(:,4);
    J = Y(:,4);
    if (length(I)==1)
        jdx = [J==1;J==2;J==3];
        Gxy = (X(I)-Y(:,I)) .* (X(J)'-Y(jdx)) ./ (Rxy.^5) .* ...
            sum((X(:,1:3)-Y(:,1:3)).*Ynrm,2);
    elseif (length(J)==1)
        idx = [I==1;I==2;I==3];
        Gxy = (X(idx)-Y(I)') .* (X(:,J)-Y(J)) ./ (Rxy.^5) .* ...
            sum((X(:,1:3)-Y(:,1:3)).*Ynrm,2);
    else
        error('oprGreenKernel.m : unavailable case')
    end
    Gxy(Rxy<1e-12) = 0;

% unknown      
else
    error('oprGreenKernel.m : unknown green kernel')
end
end
