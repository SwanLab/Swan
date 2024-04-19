%**************************************************************************
% Compliance Topology Optimization with Volume Constraint in Plane Stress,
% Steady-State Heat Conduction and Kirchhoff Plates considering multiple
% materials.
%**************************************************************************
%
% DESCRIPTION
% Computes the topological derivative and use it together with a level-set 
% domain representation method in the multi-material topology optimization
% design context.
%
% HISTORY
% S. Amstutz       06/2009: code implementation.
% A.A. Novotny     06/2009: code updating.
% D.E. Campeão     12/2010: code updating.
% J-M.C. Farias
% A.A. Novotny
% A.A. Romero      11/2018: multimaterial extension.
%**************************************************************************
% MAIN REFERENCE
% Onco, A. R., & Giusti, S. M. (2020). 
% A robust topological derivative-based multi-material optimization
% approach: Optimality condition and computational algorithm. Computer
% Methods in Applied Mechanics and Engineering, 366, 113044.
%**************************************************************************

clear; close all; clc
% load problem data
% cd('examples\elas'), [mesh, pdecoef, matprop, params, bc, psi] = viga1x2_4m;
% cd('examples\elas'), [mesh, pdecoef, matprop, params, bc, psi] = viga1x2_2m;
    
% cd('examples\elas'), [mesh, pdecoef, matprop, params, bc, psi] = viga2x1_2m_40volfrac;
cd('examples\elas'), [mesh, pdecoef, matprop, params, bc, psi] = viga2x1_3m_40volfrac;
% cd('examples\elas'), [mesh, pdecoef, matprop, params, bc, psi] = viga2x1_2m;
% cd('examples\elas'), [mesh, pdecoef, matprop, params, bc, psi] = viga2x1_3m;
% cd('examples\elas'), [mesh, pdecoef, matprop, params, bc, psi] = viga2x1_4m;

% cd('examples\heat'), [mesh, pdecoef, matprop, params, bc, psi] = engine;
% cd('examples\heat'), [mesh, pdecoef, matprop, params, bc, psi] = star_2m;
% cd('examples\heat'), [mesh, pdecoef, matprop, params, bc, psi] = star_3m;
% cd('examples\heat'), [mesh, pdecoef, matprop, params, bc, psi] = star_4m;

% cd('examples\plate'), [mesh, pdecoef, matprop, params, bc, psi] = cog_3m;
% cd('examples\plate'), [mesh, pdecoef, matprop, params, bc, psi] = cog_4m;
% cd('examples\plate'), [mesh, pdecoef, matprop, params, bc, psi] = symmetry_cog_2m_mod;
% cd('examples\plate'), [mesh, pdecoef, matprop, params, bc, psi] = symmetry_cog_2m;
% cd('examples\plate'), [mesh, pdecoef, matprop, params, bc, psi] = symmetry_cog_3m;
% cd('examples\plate'), [mesh, pdecoef, matprop, params, bc, psi] = symmetry_cog_4m;
% cd('examples\plate'), [mesh, pdecoef, matprop, params, bc, psi] = momento_2m;
% cd('examples\plate'), [mesh, pdecoef, matprop, params, bc, psi] = claDKT2;
% cd('examples\plate'), [mesh, pdecoef, matprop, params, bc, psi] = claDKT2_2m;
% cd('examples\plate'), [mesh, pdecoef, matprop, params, bc, psi] = claDKT2_3m;
% cd('examples\plate'), [mesh, pdecoef, matprop, params, bc, psi] = claDKT2_4m;

cd ..\..
filename = 'viga2x1_2m';
%Configure default settings for plot functions
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 
set(groot,'defaultAxesFontSize',16)

% mesh and geometry  parameter
p = mesh.p; e = mesh.e; t = mesh.t; g = mesh.g;
np = max(size(p)); nt = max(size(t));

% topology optimization parameters
stop = params.stop;  kmin = params.kmin; penalization = params.penalization; 
volfrac = params.volfrac; auglag = params.auglag; voleps = params.voleps;
nmat = params.nmat;
% freeze nodes
phold = [];
for i = 1:size(mesh.ghold,1)
    phold = cat(1,phold,unique(t(1:3,t(4,:)==mesh.ghold(i))));
end
pfree = setdiff(1:size(p,2),phold); % free nodes

figure(10); clf;
pdeplot(p,e,t,'xydata',t(4,:),'xystyle','flat','colormap','gray',...
              'xygrid','off','colorbar','off','title','Geometry'); 
axis image; axis off;

figure(1); clf;
set(1,'WindowStyle','docked');
pdemesh(p,e,t); axis image; axis off;

% shape functional and volume associated to the hold-all domain
psi_holdall = ones(length(p),length(matprop.E)-1); psi_holdall(:,1) = -1;
[U,F,vol_holdall] = solve(psi_holdall, mesh, matprop, pdecoef, bc); 
params.energy0 = 0.5 * dot(F,U); % initial energy 
max_vol= vol_holdall(1); %volinit = volinit.*gamma; %calculate maximum vol for each fase
params.max_vol = max_vol; % store maximum volume for each fase
voltarget = max_vol.*volfrac;

% stop criterion over the volume constraint 
volstop = voleps.*max_vol;
if penalization == 1
    volstop = max_vol;
end

[~,unitM,~] = assema(p,t,0,1,0); % mass matrix of unity density
psi = psi./normL2( unitM,psi ); % level-set function nomalization 

% plot hold-all domain
% figure(20); clf;
% set(2,'WindowStyle','docked');
% pdeplot(p,e,t,'xydata',ones(size(p,2),1),'xystyle','flat','colormap','gray',...
%               'xygrid','off','colorbar','off'); axis image; axis off; 


k = 1; iter = 0;
gsf = []; gEpot = []; gth = []; gvol = [];
option = 'null';

e_h = params.e_h;
Eps = e_h*mesh.h;
while not(strcmp(option,'s'))

    %-----------Initial values-----------%
            
    [ fi,tfi ] = charfunc( p,t,psi);
    it1 = t(1,:); it2 = t(2,:); it3 = t(3,:);
%     [K_s,M_s,F_s] = assema(p,t,Eps^2,1,tfi(1,:)) ; %before scaling
    [K_s,M_s,~] = assema(p*params.scale.L/Eps,t,1,1,0) ; %after scaling
    for jj = 1:params.nmat-1
        [~,~,F_s(:,jj)] = assema(p*params.scale.L/Eps,t,1,1,tfi(jj,:)) ; %after scaling
    end
    L = decomposition(K_s + M_s,'chol');
    tchi_s = L\F_s;
    
    [U,F,volume] = solve(psi, mesh, matprop, pdecoef, bc);
    [sf,Epot] = shfunc(F,U,volume,params);
%     sf = sf + params.beta*(2/Eps)*(1-tchi_s)'*F_s;%before scaling

    for jj = 1 : params.nmat-1
        sf = sf + params.beta*(1-tchi_s(:,jj))'*F_s(:,jj); %after scaling 
    end

    
    %%% Configure plots %%% 
    facecolor = [186 212 244]/255;
    linecolor = [0 62 102]/255;
    
    fig2 = figure(2); clf; set(2,'WindowStyle','docked'); 
    fig2.Name = 'topology';
    ploteo2 = multimat_plot( p,t,fi );
    
    fig3 = figure(3); clf; set(3,'WindowStyle','docked'); fig3.Name = 'shapefunc'; ax3 = axes();
    ploteo3{1} = plot(ax3,iter,sf,'k-','LineWidth',1.2); yyaxis left; ylabel('Shape Function'),
    yyaxis right; ploteo3{2} = plot(ax3,iter,Epot,'-','LineWidth',1.2); ylabel('Epot'),grid minor,
    
    fig4 = figure(4); clf; set(4,'WindowStyle','docked'); fig4.Name = 'thetaangle'; ax4 = axes();
    ploteo4{1} = plot(ax4,iter,0*180/pi,'-','MarkerFaceColor',facecolor,'Color',linecolor,'LineWidth',1.2); 
    hold(ax4,'on'); ploteo4{2} = plot([0 iter],[stop stop]*180/pi, '--k');
    grid minor, title('Theta Angle'); 
    
    fig5 = figure(5); clf; set(5,'WindowStyle','docked'); fig5.Name = 'volume'; ax5 = axes();
     grid minor, title('Volume Fraction');
    for index = 1:nmat-1
        ploteo5{index} = plot(ax5,iter,volume(index)/max_vol,'-','LineWidth',1.2); hold on
        leyenda{index} = sprintf('$V_{f%d}$',index);
    end
        
%     FRAMES(iter) = getframe(gcf);

    fig5 = figure(5); clf; set(5,'WindowStyle','docked'); fig5.Name = 'volume'; ax5 = axes();
    title('Volume Fraction');
    for index = 1:nmat-1
        ploteo5{index} = plot(ax5,iter,volume(index)/max_vol,'-','LineWidth',1.2); hold on
        leyenda{index} = sprintf('$V_{f%d}$',index);
    end
    grid minor

    %%%%%%%%%%%%%%%%%%%%%%%
    
    % compute function g: g = -DT (bulk) and g = DT (inclusion)
%     dtP = (2/Eps)*(1-2*tchi_s); %before scaling
    dtP = (1-2*tchi_s); %after scaling
    dtJ = topder_scaled(mesh,U,volume,matprop,psi,params,pdecoef);
    dt = dtJ + params.beta*dtP;
    dt = dt/normL2(unitM,dt);  % g function normalization
    
    dt(phold,:) = psi(phold,:); % freeze dt function
    
    for i = 1 : params.nmat-1
        cosin(i) = psi(:,i)'*unitM*dt(:,i);
    end
    cosin = max(min(sum(cosin),1.0),-1.0);
    theta = max(real(acos(cosin)),1.0e-4);
    
    difvol = volume(1:end-1)-voltarget;
    ic = (abs(difvol) > volstop); %index control
    while  (any(ic) || theta > stop || e_h > 4) && (k/2 > kmin)
                
        iter = iter + 1;
        [ fi,tfi ] = charfunc( p,t,psi);
              
        [K_s,M_s,~] = assema(p*params.scale.L/Eps,t,1,1,0) ; %after scaling
        L = decomposition(K_s + M_s,'chol');
        tchi_s = L\F_s;
        
        [U,F,volume] = solve(psi, mesh, matprop, pdecoef, bc);
        [sf,Epot] = shfunc(F,U,volume,params);
        
        for jj = 1 : params.nmat-1
            sf = sf + params.beta*(1-tchi_s(:,jj))'*F_s(:,jj); %after scaling
        end
        
        % compute function g: g = -DT (bulk) and g = DT (inclusion)
        %     dtP = (2/Eps)*(1-2*tchi_s); %before scaling
        dtP = (1-2*tchi_s); %after scaling
        dtJ = topder_scaled(mesh,U,volume,matprop,psi,params,pdecoef);
        dt = dtJ + params.beta*dtP;
        dt = dt/normL2(unitM,dt);  % g function normalization
        
        dt(phold,:) = psi(phold,:); % freeze dt function
        
        for i = 1 : length(matprop.E)-1
            cosin(i) = psi(:,i)'*unitM*dt(:,i);
        end
        cosin = max(min(sum(cosin),1.0),-1.0);
        theta = max(real(acos(cosin)),1.0e-4);
        
        % performs a line-search
        sfold = sf; psiold = psi; sf = sf + 1; k = min(1-(1e-5),1.5*k); % trick
%         
        while and((sf > sfold) , k>kmin)
%             (sf - sfold > eps(sf))
            % update level-set function
            psi(pfree,:)  = (sin((1-k)*theta)*psiold(pfree,:)...
                +  sin(k*theta)*dt(pfree,:))./sin(theta);
            [U,F,volume] = solve(psi,mesh,matprop,pdecoef,bc);
            [ ~,tfi ] = charfunc( p,t,psi);
            
            for jj = 1:params.nmat-1
                [~,~,F_s(:,jj)] = assema(p*params.scale.L/Eps,t,1,1,tfi(jj,:)) ; %after scaling
            end
            tchi_s = L\F_s;

            [sf,Epot] = shfunc(F,U,volume,params);
            for jj = 1 : params.nmat-1
                sf = sf + params.beta*(1-tchi_s(:,jj))'*F_s(:,jj); %after scaling
            end
            k = k / 2;
        end
        psi = psi/normL2(unitM,psi);
        k = k * 2;
        
        [fi,~] = charfunc( p,t,psi); %calculate characteristic function for plot
        
        gsf  = [gsf,sf];
        gEpot = [gEpot,Epot];
        gth  = [gth,theta];
        gvol = [gvol,volume(1:end-1)'*100/max_vol];
        
        disp(['iter    = ', num2str(iter)]);
        disp(['volume  = ', num2str(volume(1:end-1)),' => ', sprintf(repmat(' %3.2f%%',1,length(volume)-1),volume(1:end-1)*100/max_vol)]);
        disp(['sf      = ', num2str(sf)]);
        disp(['k       = ', num2str(k)]);
        disp(['theta   = ', num2str(theta*180/pi)]);
        disp(['penalty = ', num2str(params.penalty)]);
        disp(['e_h = ', num2str(e_h)]);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        
        facecolor = [186 212 244]/255;
        linecolor = [0 62 102]/255;
        %%% Update topology plot %%%
        drawnow
        for i = 1:nmat
            FI = [fi(it1,i)'; fi(it2,i)'; fi(it3,i)'];
            ploteo2(i).Vertices(:,3) = FI(:);
        end
%         FRAMES(iter) = getframe(fig2);
        %%% Update topology plot %%%
        
%         FRAMES(iter) = getframe(gcf);
        ploteo3{1}.XData = 1:iter; ploteo3{1}.YData = gsf;
        ploteo3{2}.XData = 1:iter; ploteo3{2}.YData = gEpot; ax3.XLim = [0 iter];
        ploteo4{1}.XData = 1:iter; ploteo4{1}.YData = gth*180/pi; ax4.XLim = [0 iter];
        ploteo4{2}.XData = [0 iter];
        for index = 1:nmat-1
            ploteo5{index}.XData = 1:iter;
            ploteo5{index}.YData = gvol(index,:);
            leyenda{index} = sprintf('$V_{f%d}$',index);
        end
        ax5.XLim = [0 iter];
        ax5.YLim = [0 105];            
        for index = 1:nmat-1
                ploteo5{index}.XData = 1:iter; 
                ploteo5{index}.YData = gvol(index,:);
                leyenda{index} = sprintf('$V_{f%d}$',index);
        end
         ax5.XLim = [0 iter];
         ax5.YLim = [0 105];
%         if params.penalization == 3
%             ploteo5{2}.XData = [0 iter]; ploteo5{3}.XData = [0 iter];
%         end



        %Update augmented lagrangian parameters
        difvol = volume(1:end-1)-voltarget;
        ic = (abs(difvol) > volstop); %index control
        
        if any(ic==1) 
            if penalization == 2 % increase the penalization parameter
                params.penalty = 2.0 * params.penalty;
                k = 1;
            elseif penalization == 3 % increase the lagrangian multiplier
                coef = volume(1:end-1)./ voltarget; coef = coef(ic); tau = auglag;
                penalty = params.penalty(ic);
%                 penalty = penalty + (tau(ic)./auglag(ic)) .* (max(0,penalty - auglag(ic).*(1.0-coef))-penalty);
                penalty = penalty + (tau(ic)./auglag(ic)) .* ((penalty - auglag(ic).*(1.0-coef))-penalty); %testing equality constraints
                params.penalty(ic) = penalty;
                
                k = 1;
            end
        end
        
        if ~all(ic) && theta < stop && e_h > 4
            e_h = e_h/2;
            Eps = e_h*mesh.h;
        end
    end
        
    option = 'null';
    while and(not(strcmp(option,'r')), and(not(strcmp(option,'s')), not(strcmp(option,'c'))))
        option = input('\n -> type "r" to remesh or "s" to stop : ', 's');
    end
    if (option == 'r')
        for i = 1:(strcmp(mesh.remesh,'longest')*2 + strcmp(mesh.remesh,'regular'))
            [p,e,t,psi] = refinemesh(g,p,e,t,psi,mesh.remesh);
            mesh.p = p; mesh.e = e; mesh.t = t;
        end
        
        % update mesh parameters
        figure(1); clf;
        set(1,'WindowStyle','docked');
        pdemesh(p,e,t); axis image; axis off;
        [~,unitM,~] = assema(p,t,0,1,0);
        area = pdetrg(p,t); mesh.area = area;
        k = 1;
        
        %Update bc
        ldir = bc.Neu.ldir;
        ldof = bc.Neu.ldof;
        lval = bc.Neu.lval;
        pbc = bouncon(p,e,t,ldir,ldof,lval);
        pbc = unique(pbc,'rows'); %just for symetry_cog_2 eliminate otherwise
        bc.pNeu = pbc;
        
        ldir = bc.Dir.ldir;
        ldof = bc.Dir.ldof;
        lval = bc.Dir.lval;
        pbc = bouncon(p,e,t,ldir,ldof,lval);
        pbc = unique(pbc,'rows'); %just for symetry_cog_2 eliminate otherwise
        bc.pDir = pbc;
        
        % freeze groups
        phold = [];
        for i = 1:size(mesh.ghold,1)
            phold = cat(1,phold,unique(t(1:3,t(4,:)==mesh.ghold(i)))); %nodes to block
        end
        pfree = setdiff(1:size(p,2),phold);
        option = 'null';
        
        % shape functional and volume associated to the hold-all domain
        psi_holdall = ones(length(p),length(matprop.E)-1); psi_holdall(:,1) = -1;
        [U,F,vol_holdall] = solve(psi_holdall, mesh, matprop, pdecoef, bc);
        params.energy0 = 0.5 * dot(F,U); % initial energy
        max_vol = vol_holdall(1);
        params.max_vol = max_vol; % store maximum volume for each fase
        voltarget = max_vol.*volfrac;
        
        % stop criterion over the volume constraint
        volstop = voleps.*max_vol;
        if penalization == 1
            volstop = max_vol;
        end
    end    
end

fig = input(' -> export figures? "y" or "n"        : ', 's');
if strcmp(fig,'y')
    cd('results')
    save(strcat(filename,'.mat'))
        print (1, '-dtiff',  strcat(filename,'_meshfem'));
        print (2, '-dtiff', strcat(filename,'_iniguess'));
%         print (3, '-djpeg', strcat(filename,'_topology'));
        saveas(3,strcat(filename,'_topology','.fig'))
        print (4, '-dtiff',  strcat(filename,'_shapefunc'));
        print (5, '-dtiff',  strcat(filename,'_thetaangle'));
        print (6, '-dtiff',  strcat(filename,'_volume'));
    cd ..
else
    return;
end