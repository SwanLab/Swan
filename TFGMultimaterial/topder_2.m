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

function dt = topder_2(mesh,U,volume,matprop,psi,params,pdecoef)

    p = mesh.p; t = mesh.t;
    penalization = params.penalization; penalty = params.penalty;
    volfrac = params.volfrac; auglag = params.auglag; volinit = params.volinit; 
    energy0 = params.energy0; voltarget = volinit*volfrac;
    %% topological derivative of the compliance shape functional
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
        
        %% topological derivative at the bulk fase
        
        for i = 1:length(matprop.k)
            gamma  = matprop.k./matprop.k(i);
            dt(:,:,i) = tk'.*((1-gamma)./(1+gamma)).*(s(1,:).*s(1,:) + s(2,:).*s(2,:))';
        end
%         dt = permute(dt,[2 1 3]);
        
    elseif size(pdecoef.f,1)==2 % vectorial second order pde
        TD = cell(length(matprop.E),length(matprop.E));
        %% material properties
        E1 = matprop.E(1); nu = matprop.nu;
        % Lame's coefficients
        la0 = nu*E1/((1+nu)*(1-2*nu)); mu0 = E1/(2*(1+nu));
        la0 = 2*mu0*la0/(la0+2*mu0); % plane stress
        gamma  = matprop.E./E1; 
        % element characteristic function: 1 if psi>0 and gama if psi<=0 
        [~,tfi] = charfunc( p,t,psi);        
%         tgamma = pdeintrp(p,t,fi*gamma'); %P1 projection method
        tgamma = gamma*tfi; %Mixed formulation method
        tE = E1*tgamma; beta = (1+nu)/(1-nu); alpha = (3-nu)/(1+nu);
        %% nominal stress
        [ux,uy]=pdegrad(p,t,U); % solution gradient
        e=[ux(1,:);(ux(2,:)+uy(1,:))/2;uy(2,:)]; % strain
        id=[1 0 1]';
        s=la0*id*(e(1,:)+e(3,:))+2*mu0*e; % stress
        % effective stress
        tgamma3 = [tgamma;tgamma;tgamma]; 
        s = s.*tgamma3;

        %% topological derivative
        %Topological derivatives are stored as following:
        % [dt1->1 dt1->2 dt1->3 ... dt1->n] -->sub matrix 1
        % [dt2->1 dt2->2 dt2->3 ... dt2->n] -->sub matrix 2
        % [dti->1 dti->2 dti->3 ... dti->n] -->sub matrix i
        % [dtn->1 dtn->2 dtn->3 ... dtn->n] -->sub matrix n
        
        % Diagonal is zero by definition of topological derivative
        % Upper triangular is dte and lower triangular is dti
        for i = 1:length(matprop.E)
            gamma  = matprop.E./matprop.E(i);
            coef1 = 0.5*((1-gamma)./(1+alpha*gamma))./tE';
            coef2 = coef1.*((gamma.*(alpha-2*beta)-1)./(1+beta*gamma));
            dt(:,:,i) = coef1'.*(4*(s(1,:).*s(1,:)+2*s(2,:).*s(2,:)+s(3,:).*s(3,:))) ...
                + coef2'.*((s(1,:)+s(3,:)).*(s(1,:)+s(3,:))) ;
        end
        dt = permute(dt,[2 1 3]);
        
        for i = 1:length(matprop.E)
            for j = 1:length(matprop.E)
                gamma  = matprop.E(j)/matprop.E(i);
                coef1 = 0.5*((1-gamma)./(1+alpha*gamma))./tE;
                coef2 = coef1.*((gamma.*(alpha-2*beta)-1)./(1+beta*gamma));
                TD{i,j} = coef1.*(4*(s(1,:).*s(1,:)+2*s(2,:).*s(2,:)+s(3,:).*s(3,:))) ...
                + coef2.*((s(1,:)+s(3,:)).*(s(1,:)+s(3,:))) ;
            end
        end
                
    elseif size(pdecoef.f,1) == 3 % scalar fourth order pde
        TD = cell(length(matprop.E),length(matprop.E));
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
        for i = 1:length(matprop.E)
            gamma  = matprop.E./matprop.E(i);
            coef0 = 6.0*(1.0 - gamma)./(tE'*h3);
            coef1 = 4.0*(alpha*beta)./(1.0 + beta*gamma);
            coef2 = 1.0./(1.0 + alpha*gamma) - 2.0*(alpha*beta)./(1.0 + beta*gamma);
            
            dt(:,:,i) = coef0.*(coef1.*aa' + coef2.*bb');
        end
%         dt = permute(dt,[2 1 3]);
        for i = 1:length(matprop.E)
            for j = 1:length(matprop.E)
                gamma  = matprop.E(j)/matprop.E(i);
                coef0 = 6.0*(1.0 - gamma)./(tE'*h3);
                coef1 = 4.0*(alpha*beta)./(1.0 + beta*gamma);
                coef2 = 1.0./(1.0 + alpha*gamma) - 2.0*(alpha*beta)./(1.0 + beta*gamma);
                TD{i,j} = coef0.*(coef1.*aa' + coef2.*bb');
            end
        end
    end
    %% include the volume constraint sensitivity 
        %Topological derivatives are stored in one matrix as following:
        % [dt1->1 dt1->2 dt1->3 ... dt1->n] 
        % [dt2->1 dt2->2 dt2->3 ... dt2->n] 
        % [dti->1 dti->2 dti->3 ... dti->n] 
        % [dtn->1 dtn->2 dtn->3 ... dtn->n] 
        
        % Diagonal is zero by definition of topological derivative
        % Upper triangular is dte and lower triangular is dti
        gamma  = matprop.E./matprop.E(1);
        for i = 1:length(matprop.E)
            for j = 1:length(matprop.E)
                dtvol(i,j) = -(gamma(i)-gamma(j)); %calculate top-der of each vol fase
            end
        end
%     lambda = [penalty'; 1];
%     dtvol = lambda.*dtvol;
%     dtvol(end,:) = lambda'.*dtvol(end,:);         %multiply for correspondig penalty parameter
       
    dt = dt / energy0;              %normalize energy topological derivative
    
    if length(matprop.E)==2
        dtE1_0 = dt(:,2,1); dtE0_1 = dt(:,1,2);
        
        dtv1_0 = dtvol(1,2); dtv0_1 = dtvol(2,1);
        
        dtv1_0 = -1; dtv0_1 = 1;
    elseif length(matprop.E)==3
        dtE1_0 = dt(:,3,1); dtE0_1 = dt(:,1,3);
        dtE2_0 = dt(:,3,2); dtE0_2 = dt(:,2,3);
        dtE1_2 = dt(:,2,1); dtE2_1 = dt(:,1,2);
        
        dtv1_0 = dtvol(1,3); dtv0_1 = dtvol(3,1);
        dtv2_0 = dtvol(2,3); dtv0_2 = dtvol(3,2);
        dtv1_2 = dtvol(1,2); dtv2_1 = dtvol(2,1);
        
        dtv1_0 = -1; dtv0_1 = 1;
        dtv2_0 = -1; dtv0_2 = 1;
        dtv1_2 = -1; dtv2_1 = 1;
        
    elseif length(matprop.E)==4
        dtE1_0 = dt(:,4,1); dtE0_1 = dt(:,1,4);
        dtE2_0 = dt(:,4,2); dtE0_2 = dt(:,2,4);
        dtE3_0 = dt(:,4,3); dtE0_3 = dt(:,3,4);
        dtE1_2 = dt(:,2,1); dtE2_1 = dt(:,1,2);
        dtE1_3 = dt(:,3,1); dtE3_1 = dt(:,1,3);
        dtE2_3 = dt(:,3,2); dtE3_2 = dt(:,2,3);
        
        dtv1_0 = dtvol(1,4); dtv0_1 = dtvol(4,1);
        dtv2_0 = dtvol(2,4); dtv0_2 = dtvol(4,2);
        dtv3_0 = dtvol(3,4); dtv0_3 = dtvol(4,3);
        dtv1_2 = dtvol(1,2); dtv2_1 = dtvol(2,1);
        dtv1_3 = dtvol(1,3); dtv3_1 = dtvol(3,1);
        dtv2_3 = dtvol(2,3); dtv3_2 = dtvol(3,2);
        
        dtv1_0 = -1; dtv0_1 = 1;
        dtv2_0 = -1; dtv0_2 = 1;
        dtv3_0 = -1; dtv0_3 = 1;
        dtv1_2 = -1; dtv2_1 = 1;
        dtv1_3 = -1; dtv3_1 = 1;
        dtv2_3 = -1; dtv3_2 = 1;
    end
    
    
    %Sum volume constraint
    if penalization == 1
        if length(matprop.E)==2
            dt1_0 = dtE1_0 + penalty(1)*dtv1_0/volinit; dt0_1 = dtE0_1 + penalty(1)*dtv0_1/volinit;       
        elseif length(matprop.E)==3
            dt1_0 = dtE1_0 + penalty(1)*dtv1_0/volinit; dt0_1 = dtE0_1 + penalty(1)*dtv0_1/volinit;
            dt2_0 = dtE2_0 + penalty(2)*dtv2_0/volinit; dt0_2 = dtE0_2 + penalty(2)*dtv0_2/volinit;
            dt1_2 = dtE1_2 + (penalty(1)/volinit - penalty(2)/volinit)*dtv1_2;
            dt2_1 = dtE2_1 + (penalty(1)/volinit - penalty(2)/volinit)*dtv2_1;
        elseif length(matprop.E)==4
            dt1_0 = dtE1_0 + penalty(1)*dtv1_0/volinit; dt0_1 = dtE0_1 + penalty(1)*dtv0_1/volinit;
            dt2_0 = dtE2_0 + penalty(2)*dtv2_0/volinit; dt0_2 = dtE0_2 + penalty(2)*dtv0_2/volinit;
            dt3_0 = dtE3_0 + penalty(3)*dtv3_0/volinit; dt0_3 = dtE0_3 + penalty(3)*dtv0_3/volinit;
            dt1_2 = dtE1_2 + (penalty(1)/volinit - penalty(2)/volinit)*dtv1_2;
            dt2_1 = dtE2_1 + (penalty(1)/volinit - penalty(2)/volinit)*dtv2_1;
            dt1_3 = dtE1_3 + (penalty(1)/volinit - penalty(3)/volinit)*dtv1_3;
            dt3_1 = dtE3_1 + (penalty(1)/volinit - penalty(3)/volinit)*dtv3_1;
            dt2_3 = dtE2_3 + (penalty(2)/volinit - penalty(3)/volinit)*dtv2_3;
            dt3_2 = dtE3_2 + (penalty(2)/volinit - penalty(3)/volinit)*dtv3_2;
        end
    elseif penalization == 3
        if length(matprop.E)==2
            coef = volume(1:end-1) ./ voltarget;
            dt1_0 = dtE1_0 + max(0,penalty(1) - auglag(1)*(1.0-coef(1))) *dtv1_0/ voltarget(1);
            dt0_1 = dtE0_1 + max(0,penalty(1) - auglag(1)*(1.0-coef(1))) *dtv0_1/ voltarget(1);
            
        elseif length(matprop.E)==3
            coef = volume(1:end-1) ./ voltarget;
            dt1_0 = dtE1_0 + max(0,penalty(1) - auglag(1)*(1.0-coef(1))) *dtv1_0/ voltarget(1);
            dt0_1 = dtE0_1 + max(0,penalty(1) - auglag(1)*(1.0-coef(1))) *dtv0_1/ voltarget(1);
            dt2_0 = dtE2_0 + max(0,penalty(2) - auglag(2)*(1.0-coef(2))) *dtv2_0/ voltarget(2);
            dt0_2 = dtE0_2 + max(0,penalty(2) - auglag(2)*(1.0-coef(2))) *dtv0_2/ voltarget(2);
            dt1_2 = dtE1_2 + (max(0,penalty(1) - auglag(1)*(1.0-coef(1)))/ voltarget(1) -  max(0,penalty(2) - auglag(2)*(1.0-coef(2))) / voltarget(2))*dtv1_2;
            dt2_1 = dtE2_1 + (max(0,penalty(1) - auglag(1)*(1.0-coef(1)))/ voltarget(1) -  max(0,penalty(2) - auglag(2)*(1.0-coef(2))) / voltarget(2))*dtv2_1;
            
        elseif length(matprop.E)==4
            coef = volume(1:end-1) ./ voltarget;
            dt1_0 = dtE1_0 + max(0,penalty(1) - auglag(1)*(1.0-coef(1))) *dtv1_0/ voltarget(1);
            dt0_1 = dtE0_1 + max(0,penalty(1) - auglag(1)*(1.0-coef(1))) *dtv0_1/ voltarget(1);
            dt2_0 = dtE2_0 + max(0,penalty(2) - auglag(2)*(1.0-coef(2))) *dtv2_0/ voltarget(2);
            dt0_2 = dtE0_2 + max(0,penalty(2) - auglag(2)*(1.0-coef(2))) *dtv0_2/ voltarget(2);
            dt3_0 = dtE3_0 + max(0,penalty(3) - auglag(3)*(1.0-coef(3))) *dtv3_0/ voltarget(3);
            dt0_3 = dtE0_3 + max(0,penalty(3) - auglag(3)*(1.0-coef(3))) *dtv0_3/ voltarget(3);
            dt1_2 = dtE1_2 + (max(0,penalty(1) - auglag(1)*(1.0-coef(1)))/ voltarget(1) -  max(0,penalty(2) - auglag(2)*(1.0-coef(2))) / voltarget(2))*dtv1_2;
            dt2_1 = dtE2_1 + (max(0,penalty(1) - auglag(1)*(1.0-coef(1)))/ voltarget(1) -  max(0,penalty(2) - auglag(2)*(1.0-coef(2))) / voltarget(2))*dtv2_1;
            dt1_3 = dtE1_3 + (max(0,penalty(1) - auglag(1)*(1.0-coef(1)))/ voltarget(1) -  max(0,penalty(3) - auglag(3)*(1.0-coef(3))) / voltarget(3))*dtv1_3;
            dt3_1 = dtE3_1 + (max(0,penalty(1) - auglag(1)*(1.0-coef(1)))/ voltarget(1) -  max(0,penalty(3) - auglag(3)*(1.0-coef(3))) / voltarget(3))*dtv3_1;
            dt2_3 = dtE2_3 + (max(0,penalty(2) - auglag(2)*(1.0-coef(2)))/ voltarget(2) -  max(0,penalty(3) - auglag(3)*(1.0-coef(3))) / voltarget(3))*dtv2_3;
            dt3_2 = dtE3_2 + (max(0,penalty(2) - auglag(2)*(1.0-coef(2)))/ voltarget(2) -  max(0,penalty(3) - auglag(3)*(1.0-coef(3))) / voltarget(3))*dtv3_2;            
        end
    end
        
    %% smoothing of the topological derivative
    
    if length(matprop.E)==2
        [~,tfi] = charfunc( p,t,psi); tfi = tfi';   % get characteristic function - Mixed formulation method
        %     fi = (pdeintrp(p,t,fi)).'; % interpolation at gauss point - P1 projection method
       
        % compute g for each level set g: g = -DT (bulk) and g = DT (inclusion)

        dt = [];
        dt = - tfi(:,1).*dt1_0 + tfi(:,end).*dt0_1;

    elseif length(matprop.E)==3
        [~,tfi] = charfunc( p,t,psi); tfi = tfi';   % get characteristic function - Mixed formulation method
        [tXi2,~] = integ_exact(t,p,psi(:,2)); chi2 = (1 - tXi2)'; %- Mixed formulation method
        %     fi = (pdeintrp(p,t,fi)).'; % interpolation at gauss point - P1 projection method
        %     chi2 = (pdeintrp(p,t,(psi(:,2)<0))).'; %- P1 projection method
        
        % compute g for each level set g: g = -DT (bulk) and g = DT (inclusion)
        
        dt = [];
        dt(:,1) = - tfi(:,1).*dt1_0 - tfi(:,2).*dt2_0 + tfi(:,3).*( (1-chi2).*dt0_1 + chi2.*dt0_2 );
        dt(:,2) = tfi(:,1).*dt1_2 - tfi(:,2).*dt2_1;
        
    elseif length(matprop.E)==4
        [~,tfi] = charfunc( p,t,psi); tfi = tfi';   % get characteristic function - Mixed formulation method
        [tXi2,~] = integ_exact(t,p,psi(:,2)); chi2 = (1 - tXi2)'; %- Mixed formulation method
        [tXi3,~] = integ_exact(t,p,psi(:,3)); chi3 = (1 - tXi3)'; %- Mixed formulation method
        %     fi = (pdeintrp(p,t,fi)).'; % interpolation at gauss point - P1 projection method
        %     chi2 = (pdeintrp(p,t,(psi(:,2)<0))).'; %- P1 projection method
        
        % compute g for each level set g: g = -DT (bulk) and g = DT (inclusion)

        dt = [];
        dt(:,1) = - tfi(:,1).*dt1_0 - tfi(:,2).*dt2_0 - tfi(:,3).*dt3_0+ tfi(:,4).*( (1-chi2).*dt0_1 + (1-chi3).*chi2.*dt0_2 + chi2.*chi3.*dt0_3 );
        dt(:,2) = - tfi(:,2).*dt2_1 - tfi(:,3).*dt3_1 + tfi(:,1).*( (1-chi3).*dt1_2 + chi3.*dt1_3 );
        dt(:,3) = tfi(:,2).*dt2_3 - tfi(:,3).*dt3_2;

    end
    
         dt = pdeprtni(p,t,dt.'); %P1-projection method
%    [ dt ] = mixedpdeprtni( p,t,dt' ); %mixed formulation aproach
    
end