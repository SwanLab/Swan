%**************************************************************************
% Topological Derivative for the Compliance with Volume Constraint  
%**************************************************************************
% DESCRIPTION
% Computes the topological derivative
% 
% INPUT
% mesh:    pdetool mesh struct
% U:       pdetool solution
% volume:  current volume
% matprop: material properties struct
% psi:     level-set function
% params:  topology optimization parameters struct
% pdecoef: pdetool coeficiets struct
%
% OUTPUT
% dt: topological derivative -> dt = -DT (bulk) and dt = DT (inclusion)
%
% HISTORY
% S. Amstutz     06/2009: code implementation.
% A.A. Novotny   06/2009: code updating.
% D.E. Campeão   12/2010: code updating.
% J-M.C. Farias
% A.A. Novotny
%**************************************************************************

function dt = topder(mesh,U,volume,matprop,psi,params,pdecoef,designVariable,m, TOParams)

    p = mesh.p; t = mesh.t;
    penalization = params.penalization; penalty = params.penalty;
    volfrac = params.volfrac; auglag = params.auglag; max_vol = TOParams.max_vol; 
    energy0 = TOParams.energy0; voltarget = max_vol*volfrac;
    nmat = 4; % changed 
    %% topological derivative of the compliance shape functional
    % Topological derivatives are stored as following:
    % TD{i,j} ------> dti->j (perturbatation made of material "j" inside the material domain "i")
    % Diagonal is zero by definition of topological derivative
    % Upper triangular is dte and lower triangular is dti
    
    if size(pdecoef.f,1)==1 % scalar second order pde 
        TD = cell(length(matprop.k),length(matprop.k));
        %% material properties
        k1 = matprop.k(1);
        gamma = matprop.k/matprop.k(1);
        % element characteristic function: 1 if psi>0 and gama if psi<=0
        [~,tfi] = charfunc( p,t,psi); 
        tgamma = gamma*tfi; %Mixed formulation method
        tk = k1*tgamma;    
        [ux,uy] = pdegrad(p,t,U); % solution gradient
        s = [ux(1,:);uy(1,:)]; % nominal heat flux
        
        %% topological derivative    
        
        for i = 1:nmat
            for j = 1:nmat
                gamma  = matprop.k(j)./matprop.k(i);
                TD{i,j} = tk.*((1-gamma)./(1+gamma)).*(s(1,:).*s(1,:) + s(2,:).*s(2,:));
            end
        end
        
    elseif size(pdecoef.f,1)==2 % vectorial second order pde
        TD = cell(nmat,nmat);
        %% material properties
        E1 = matprop.A.young; nu = matprop.A.nu;
        % Lame's coefficients
        la0 = nu.*E1./((1+nu).*(1-2.*nu)); mu0 = E1./(2.*(1+nu)); % plane strain
        la0 = 2.*mu0.*la0./(la0+2.*mu0); % plane stress
        gamma(1)  = matprop.A.young/E1;
        gamma(2)  = matprop.B.young/E1; 
        gamma(3)  = matprop.C.young/E1; 
        gamma(4)  = matprop.D.young/E1; 
        % element characteristic function: 1 if psi>0 and gama if psi<=0 
        cParams.psi = psi;
        cParams.p = p;
        cParams.t = t; 
        cParams.designVariable = designVariable;
        cParams.m = m;
        charfun = CharacteristicFunctionComputer(cParams); % s'ha de construir la classe - charfunc!!
        [~,tfi] = charfun.computeFiandTfi();
%         tgamma = pdeintrp(p,t,fi*gamma'); %P1 projection method
        tgamma = gamma*tfi; %Mixed formulation method
        tE = E1*tgamma; beta = (1+nu)/(1-nu); alpha = (3-nu)/(1+nu);
        E = [matprop.A.young matprop.B.young matprop.C.young matprop.D.young];
        %% nominal stress
        [ux,uy]=pdegrad(p,t,U); % solution gradient
        e=[ux(1,:);(ux(2,:)+uy(1,:))/2;uy(2,:)]; % strain
        id=[1 0 1]';
        s=la0*id*(e(1,:)+e(3,:))+2*mu0*e; % stress
        % effective stress
        tgamma3 = [tgamma;tgamma;tgamma]; 
        s = s.*tgamma3;

        %% topological derivative
        
        for i = 1:nmat
            for j = 1:nmat
                gamma  = E(j)/E(i);
                coef1 = 0.5*((1-gamma)./(1+alpha*gamma))./tE;
                coef2 = coef1.*((gamma.*(alpha-2*beta)-1)./(1+beta*gamma));
                TD{i,j} = coef1.*(4*(s(1,:).*s(1,:)+2*s(2,:).*s(2,:)+s(3,:).*s(3,:))) ...
                + coef2.*((s(1,:)+s(3,:)).*(s(1,:)+s(3,:))) ;
            end
        end
                
    elseif size(pdecoef.f,1) == 3 % scalar fourth order pde
        TD = cell(nmat,nmat);
        %% material properties
        h0 = matprop.h0; h3 = h0*h0*h0;
        E1 = matprop.E(1); nu = matprop.nu;
        gamma = matprop.E./matprop.E(1); % contrast 
        % element characteristic function: 1 if psi>0 and gama if psi<=0
        [~,tfi] = charfunc( p,t,psi);  
        tgamma = gamma*tfi; % Mixed formulation method 
        tE = E1*tgamma; alpha = (1+nu)/(1-nu); beta = (1-nu)/(3+nu);
        %% nominal stress
        s = getstress(U,mesh,matprop,pdecoef);
        % effective stress
        tgamma3 = [tgamma;tgamma;tgamma]; 
        s = s.*tgamma3;
        aa = s(1,:).*s(1,:) + 2*s(2,:).*s(2,:) + s(3,:).*s(3,:); % s . s
        bb = (s(1,:) + s(3,:)).*(s(1,:) + s(3,:)); % (tr(s))^2

        %% topological derivative         
        for i = 1:nmat
            for j = 1:nmat
                gamma  = matprop.E(j)/matprop.E(i);
                coef0 = 6.0*(1.0 - gamma)./(tE*h3);
                coef1 = 4.0*(alpha*beta)./(1.0 + beta*gamma);
                coef2 = 1.0./(1.0 + alpha*gamma) - 2.0*(alpha*beta)./(1.0 + beta*gamma);
                TD{i,j} = coef0.*(coef1.*aa + coef2.*bb);
            end
        end
    end
    %% include the volume constraint sensitivity          
    coef = volume(1:end-1) ./ voltarget;
    for i = 1:nmat
        for j = 1:nmat
            if penalization == 1
                if i==j
                    TD{i,j} = 0;
                elseif j == nmat
                    TD{i,j} = TD{i,j}/energy0 - penalty(i)/max_vol;
                elseif i == nmat
                    TD{i,j} = TD{i,j}/energy0 + penalty(j)/max_vol;
                else
                    TD{i,j} = TD{i,j}/energy0 + (penalty(j) - penalty(i))/max_vol;
                end
            elseif penalization == 3
                if i==j
                    TD{i,j} = 0;
                elseif j == nmat
                    TD{i,j} = TD{i,j}/energy0 - ( (penalty(i) + auglag(i)*(coef(i)-1))/voltarget(i) );
                elseif i == nmat
                    TD{i,j} = TD{i,j}/energy0 + ( (penalty(j) + auglag(j)*(coef(j)-1))/voltarget(j) );
                else
                    TD{i,j} = TD{i,j}/energy0 + ( (penalty(j) + auglag(j)*(coef(j)-1))/voltarget(j) )...
                                              - ( (penalty(i) + auglag(i)*(coef(i)-1))/voltarget(i) );
                end
            end
        end
    end
    %% smoothing of the topological derivative
    
    if nmat==2
        [~,tfi] = charfunc( p,t,psi); % get characteristic function - Mixed formulation method
        %     fi = (pdeintrp(p,t,fi)).'; % interpolation at gauss point - P1 projection method
       
        % compute g for each level set g: g = -DT (bulk) and g = DT (inclusion)

        dt = [];
        dt = - tfi(1,:).*TD{1,end} + tfi(end,:).*TD{end,1};

    elseif nmat==3
        [~,tfi] = charfunc( p,t,psi);  % get characteristic function - Mixed formulation method
        [tXi2,~] = integ_exact(t,p,psi(:,2)); chi2 = (1 - tXi2); %- Mixed formulation method
        %     fi = (pdeintrp(p,t,fi)).'; % interpolation at gauss point - P1 projection method
        %     chi2 = (pdeintrp(p,t,(psi(:,2)<0))).'; %- P1 projection method
        
        % compute g for each level set g: g = -DT (bulk) and g = DT (inclusion)
        
        dt = [];
        dt(1,:) = - tfi(1,:).*TD{1,end} - tfi(2,:).*TD{2,end} + tfi(3,:).*( (1-chi2).*TD{end,1} + chi2.*TD{end,2} );
        dt(2,:) = tfi(1,:).*TD{1,2} - tfi(2,:).*TD{2,1};
        
    elseif nmat==4
        charfun = CharacteristicFunctionComputer(cParams); 
        [~,tfi] = charfun.computeFiandTfi();
        [tXi2,~] = integ_exact(t,p,psi(:,2)); chi2 = (1 - tXi2); %- Mixed formulation method
        [tXi3,~] = integ_exact(t,p,psi(:,3)); chi3 = (1 - tXi3); %- Mixed formulation method
        %     fi = (pdeintrp(p,t,fi)).'; % interpolation at gauss point - P1 projection method
        %     chi2 = (pdeintrp(p,t,(psi(:,2)<0))).'; %- P1 projection method
        
        % compute g for each level set g: g = -DT (bulk) and g = DT (inclusion)

        dt = [];
        dt(1,:) = - tfi(1,:).*TD{1,end} - tfi(2,:).*TD{2,end} - tfi(3,:).*TD{3,end} ...
            + tfi(4,:).*( (1-chi2).*TD{end,1} + (1-chi3).*chi2.*TD{end,2} + chi2.*chi3.*TD{end,3} );
        dt(2,:) = - tfi(2,:).*TD{2,1} - tfi(3,:).*TD{3,1} + tfi(1,:).*( (1-chi3).*TD{1,2} + chi3.*TD{1,3} );
        dt(3,:) = tfi(2,:).*TD{2,3} - tfi(3,:).*TD{3,2};

    end
    
         dt = pdeprtni(p,t,dt); %P1-projection method
%    [ dt ] = mixedpdeprtni( p,t,dt ); %mixed formulation aproach
    
end