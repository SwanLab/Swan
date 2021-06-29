function Ms = domRegularize2D(data)
%+========================================================================+
%|                                                                        |
%|              OPENDOM - LIBRARY FOR NUMERICAL INTEGRATION               |
%|           openDom is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : domRegularize2D.m                             |
%|    #    |   VERSION    : 0.53                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 25.11.2018                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2019                                    |
%| ( === ) |   SYNOPSIS   : Finite element regularization matrix for      |
%|  `---'  |                singularities with Laplace kernel             |
%+========================================================================+

%%% INPUT ANALYSIS
if (length(data) == 4)
    X     = data{1};
    Ydom  = data{2};
    green = data{3};
    v     = data{4};
else
    Xdom  = data{1};
    Ydom  = data{2};
    u     = data{3};
    green = data{4};
    v     = data{5};
end


%%% INITIALIZATION
% Mesh data from Y
vtx  = Ydom.msh.vtx;
elt  = Ydom.msh.elt;
ctr  = Ydom.msh.ctr;
nrm  = Ydom.msh.nrm;
stp  = Ydom.msh.stp;
lgt  = Ydom.msh.ndv;
tau  = Ydom.msh.tgt;
nu   = cell2mat(Ydom.msh.nrmEdg);
Nelt = size(elt,1);

% Quadrature data from Y
[Y,Wy,elt2qud] = Ydom.qud;
Yun            = ones(1,size(elt2qud,2));

% Degrees of freedom from Y
[~,elt2dof] = v.dof;
Nbas        = size(elt2dof,2);

% Quadrature data from X
if (length(data) == 5)
    [X,Wx] = Xdom.qud;
end
Nx = size(X,1);

% Rangesearch with max(|edge|)_Y
[Ielt,Relt] = rangeSearch(X,ctr,1.1*stp(2));                  %%% DEBUG %%%
Mx          = cell(Nelt,1);


%%% RIGHT INTEGRATION WITH REGULARIZATION
for el = 1:Nelt
    % Edge data for Y
    Sel  = vtx(elt(el,:),:);
    Nel  = nrm(el,:);
    Tel  = tau(el,:);
    NUel = reshape(nu(el,:),3,2)';
    
    % Local size
    rMin = 1.1*lgt(el);                                       %%% DEBUG %%%
    
    % Quadratures points in interaction
    Iy = elt2qud(el,:);
    Ix = sort(Ielt{el}(Relt{el}<rMin))';
    
%     % Graphical representation                              %%% DEBUG %%%
%     figure(10)
%     plot(Ydom.msh.sub(el))
%     hold on
%     plot3(X(Ix,1),X(Ix,2),X(Ix,3),'*')
%     axis equal
%     hold off
%     pause()
    
    % If interactions
    if ~isempty(Ix)
        %%% CORRECTION WITH SEMI-ANALYTIC INTEGRATION
        % Analytical integration
        [logR,rlogR,gradlogR] = domSemiAnalyticInt2D(X(Ix,:),Sel,Nel,Tel);
%         logR(:) = 0; rlogR(:) = 0; gradlogR(:) = 0;                                         %%% DEBUG %%%
        
        % Vector yg-x
        Xun = ones(length(Ix),1);
        XY1 = Xun * Y(Iy,1)' - X(Ix,1) * Yun;
        XY2 = Xun * Y(Iy,2)' - X(Ix,2) * Yun;
        
        % Distance r = |yg-x|
        Rxy               = sqrt(XY1.^2 + XY2.^2);
        logRxy            = log(Rxy);
        logRxy(Rxy<1e-12) = 0;
        
        % Int_el(log(|r|)) - Sum_g log(|yg-x|)
        logR = logR - logRxy * Wy(Iy);
%         norm(logR)                                          %%% DEBUG %%%
                
        % Int_el(r*log(|r|)) - Sum_g (yg-x)*log(|yg-x|)
        rlogR(:,1) = rlogR(:,1) - (XY1 .* logRxy) * Wy(Iy);
        rlogR(:,2) = rlogR(:,2) - (XY2 .* logRxy) * Wy(Iy);
        rlogR(:,3) = 0;
%         norm(rlogR)                                         %%% DEBUG %%%
        
        % Int_el(r/|r|^2) - Sum_g (yg-x)/|yg-x|^2
        Rxym12            = 1./(Rxy.^2);
        Rxym12(Rxy<1e-12) = 0;
        gradlogR(:,1)     = gradlogR(:,1) - (XY1 .* Rxym12) * Wy(Iy);
        gradlogR(:,2)     = gradlogR(:,2) - (XY2 .* Rxym12) * Wy(Iy);
        gradlogR(:,3)     = 0;
%         norm(gradlogR)                                      %%% DEBUG %%%
        
        % Nullify V
        V = [];
        
        %%% FINITE ELEMENT P0
        if strcmp(v.typ,'P0')
            % Correction
            if strcmp(green,'[log(r)]') && strcmp(v.opr,'[psi]')
                V = logR;       
                
            elseif strcmp(green,'grady[log(r)]') && strcmp(v.opr,'n*[psi]')
                V = gradlogR * Nel';   
                
            elseif strcmp(green,'[log(r)]') && strcmp(v.opr,'n*[psi]')
                V{1} = logR .* Nel(1);
                V{2} = logR .* Nel(2);
                V{3} = logR .* Nel(3);
                
            elseif strcmp(green,'[log(r)]') && strcmp(v.opr,'nxgrad[psi]')
                V{1} = zeros(size(logR));
                V{2} = zeros(size(logR));
                V{3} = zeros(size(logR));
                
            else
                error('domRegularize2D.m : unavailable case')
            end
            
        %%% FINITE ELEMENT P1
        elseif strcmp(v.typ,'P1')
            % For each basis function
            for j = 1:Nbas
                % Next dof
                jp1 = mod(j,2) + 1;
                
                % Height from j
                hj = ((Sel(j,:)-Sel(jp1,:)) * NUel(j,:)');
                
                % Scalar product (x-yk).nuj/hj
                tmp = ( (X(Ix,1)-Sel(jp1,1))*NUel(j,1) + ...
                    (X(Ix,2)-Sel(jp1,2))*NUel(j,2) ) ./ hj;
                
                % Correction
                if strcmp(green,'[log(r)]') && strcmp(v.opr,'[psi]')
                    V(:,j) = logR.*tmp + rlogR*NUel(j,:)'/hj;  
                    
                elseif strcmp(green,'[log(r)]') && strcmp(v.opr,'n*[psi]')
                    Vx        = logR.*tmp + rlogR*NUel(j,:)'/hj;
                    V{1}(:,j) = Vx .* Nel(1);
                    V{2}(:,j) = Vx .* Nel(2);
                    V{3}(:,j) = Vx .* Nel(3);
                    
                elseif strcmp(green,'[log(r)]') && strcmp(v.opr,'nxgrad[psi]')
                    NxNUj     = cross(Nel,NUel(j,:));
                    V{1}(:,j) = NxNUj(1)/hj .* logR;
                    V{2}(:,j) = NxNUj(2)/hj .* logR;
                    V{3}(:,j) = NxNUj(3)/hj .* logR;
                    
                elseif strcmp(green,'grady[log(r)]') && strcmp(v.opr,'n*[psi]')
                    V(:,j) = tmp .* (gradlogR * Nel');    
                    
                else
                    error('domRegularize2D.m : unavailable case')
                end
            end
            
        else
            error('domRegularize2D.m : unavailable case')
        end
        
        % Matrix-Vector product
        I = repmat(Ix,1,Nbas);
        J = repmat(elt2dof(el,:),length(Ix),1);
        if iscell(V)
            Mx{el} = [I(:) J(:) V{1}(:) V{2}(:) V{3}(:)];
        else
            Mx{el} = [I(:) J(:) V(:)];
        end
    end
end


%%% LEFT INTEGRATION AND RIGHT REDUCTION
% Left integration matrix
Ndof = size(v.dof,1);
if (length(data) == 4)
    Mu = speye(Nx,Nx);
    Mw = speye(Nx,Nx);
    Ms = sparse(Nx,Ndof);
else
    Mu = u.uqm(Xdom);
    Mw = spdiags(Wx,0,length(Wx),length(Wx));
    Ms = sparse(length(u),Ndof);
end

% Regularization matrix
Mx = double(cell2mat(Mx));
if isempty(Mx)
    if iscell(Mu)
        Mx = [1 1 zeros(1,length(Mu))];
    else
        Mx = [1 1 0];
    end
end

% Left integration
if iscell(Mu)
    for i = 1:length(Mu)
        Ms = Ms + Mu{i}' * Mw * sparse(Mx(:,1),Mx(:,2),Mx(:,2+i),Nx,Ndof);
    end
else
    Ms = Mu' * Mw * sparse(Mx(:,1),Mx(:,2),Mx(:,3),Nx,Ndof);
end

% Right reduction
[~,Mv] = v.unk;
Ms     = Ms * Mv;
end
