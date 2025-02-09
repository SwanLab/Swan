function ray = rayCollision(ray,Iray)
%+========================================================================+
%|                                                                        |
%|           OPENRAY - LIBRARY FOR TRI-DIMENSIONAL RAY TRACING            |
%|           openRay is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : rayCollision.m                                |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Intersection between ray and mesh             |
%|  `---'  |                                                              |
%+========================================================================+
 
% References : 
% https://fr.wikipedia.org/wiki/Plan_%28math%C3%A9matiques%29
% https://fr.wikipedia.org/wiki/Matrice_inversible#Inversion_des_matrices_3_x_3
% https://fr.wikipedia.org/wiki/Algorithme_de_Gram-Schmidt
%
% Plan             : OA + alpha*E1 + beta*E2
% Droite           : OB + gamma*U
% Intersection     : alpha*E1 + beta*E2 = AB + gamma*U 
% Linear system : 
%  | E1x  E2x  -Ux | * | alpha | = | Bx - Ax |
%  | E1y  E2y  -Uy |   | beta  |   | By - Ay |
%  | E1z  E2z  -Uz |   | gamma |   | Bz - Az |
 
% Element origin
A = ray.msh.vtx(ray.msh.elt(:,1),:);
 
% Element tangential basis 
E1 = ray.msh.vtx(ray.msh.elt(:,2),:) - A;
E2 = ray.msh.vtx(ray.msh.elt(:,3),:) - A;
 
% Ray origin
B = ray.pos{end};
 
% Ray basis
U = ray.dir;

% Initialiation
Dray = 1e8*ones(length(ray),1);
Ielt = 1e8*ones(length(ray),1);

% For each element
for el = 1:length(Iray)
    % Ray indices
    ir = Iray{el};
    
    % For non-empty indices
    if ~isempty(ir)
        % Right hand side
        Vx = B(ir,1) - A(el,1);
        Vy = B(ir,2) - A(el,2);
        Vz = B(ir,3) - A(el,3);
        
        % Linear system
        a = E1(el,1);
        b = E2(el,1);
        c = -U(ir,1);
        d = E1(el,2);
        e = E2(el,2);
        f = -U(ir,2);
        g = E1(el,3);
        h = E2(el,3);
        i = -U(ir,3);
        
        % Determinant
        det = a.*e.*i + b.*f.*g + c.*d.*h ...
            - c.*e.*g - f.*h.*a - i.*b.*d ;
        
        % Extract singularities
        det(abs(det)<1e-8) = 1e-8;
        
        % Exact resolution with Sarrus rules
        alpha = 1./det .* ( (e.*i-f.*h).*Vx + (c.*h-b.*i).*Vy + (b.*f-c.*e).*Vz );
        beta  = 1./det .* ( (f.*g-d.*i).*Vx + (a.*i-c.*g).*Vy + (c.*d-a.*f).*Vz );
        gamma = 1./det .* ( (d.*h-e.*g).*Vx + (b.*g-a.*h).*Vy + (a.*e-b.*d).*Vz );

        % Ray inside current element
        ind = find( (gamma > 1e-8)     & ...   % ray in the right sens
            (alpha >= 0) & (alpha <= 1) & ...  % ray coordinates in [0,1]*E1
            (beta >= 0)  & (beta  <= 1) & ...  % ray coordinates in [0,1]*E2
            (alpha+beta <= 1)  & ...           % ray inside triangle
            (gamma < Dray(ir)) ) ;             % Distance is the smallest
        
        % Update data
        if ~isempty(ind)
            Dray(ir(ind)) = gamma(ind);
            Ielt(ir(ind)) = el;
        end
    end
end

% Available ray indices
I    = (Dray < 1e8);
Ielt = Ielt(I);
Dray = Dray(I);
U    = U(I,:);

% Update positions
ray.pos{1,end+1}    = B;
ray.pos{1,end}(I,:) = B(I,:) + (Dray*[1 1 1]).*U;
 
% Update element
ray.iel(:,end+1) = ray.iel(:,end);
ray.iel(I,end)   = Ielt;
 
% Orthonormal basis for each ray  (gramm schmidt)
E1 = E1(Ielt,:);
E1 = E1 ./ (sqrt(sum(E1.^2,2)) * [1 1 1]);
 
E2 = E2(Ielt,:);
E2 = E2 - (sum(E1.*E2,2)*[1 1 1]).*E1;
E2 = E2 ./ (sqrt(sum(E2.^2,2)) * [1 1 1]);
 
E3 = cross(E1,E2,2);
 
% Update direction with normal reflexions
ray.dir(I,:) = (sum(U.*E1,2)*[1 1 1]) .* E1 ...
     + (sum(U.*E2,2)*[1 1 1]) .* E2 ...
     - (sum(U.*E3,2)*[1 1 1]) .* E3 ;

 % Verify unitarity
if (norm(sqrt(sum(ray.dir(I,:).^2,2)) - 1,'inf') > 1e-12)
    error('rayCollision.m : non unitary directions')
end  
 
% Kill unavailable directions
ray.dir(~I,:) = 0 ;
end
 