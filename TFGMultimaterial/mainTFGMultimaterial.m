function psi = mainTFGMultimaterial()


clear all; close all; clc
% load problem data

% Set penalization, line search and stop criterion parameters
params = ParametersComputer();

% Generate mesh
mesh = MeshComputer();

% Generate boundary conditions
cParams.e = mesh.e;
cParams.p = mesh.p;
cParams.t = mesh.t;
bc = BoundaryConditionsComputer(cParams);


[mesh, pdecoef, matprop, params, bc, psi] = viga2x1_4m;


filename = 'viga2x1_2m';
%Configure default settings for plot functions
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 
set(groot,'defaultAxesFontSize',16)

% mesh and geometry  parameter
p = mesh.p; e = mesh.e; t = mesh.t; g = mesh.g;

% topology optimization parameters
stop = params.stop;  kmin = params.kmin; penalization = params.penalization; 
volfrac = params.volfrac; auglag = params.auglag; voleps = params.voleps;
% freeze nodes
phold = [];
for i = 1:size(mesh.ghold,1)
    phold = cat(1,phold,unique(t(1:3,t(4,:)==mesh.ghold(i))));
end
pfree = setdiff(1:size(p,2),phold); % free nodes

%figure(10); clf;
pdeplot(p,e,t,'xydata',t(4,:),'xystyle','flat','colormap','gray',...
              'xygrid','off','colorbar','off','title','Geometry'); 
axis image; axis off;

%figure(1); clf;
%set(1,'WindowStyle','docked');
pdemesh(p,e,t); axis image; axis off;

% shape functional and volume associated to the hold-all domain
psi_hold_all = ones(length(p),length(matprop.E)-1); psi_hold_all(:,1) = -1;
[U,F,vol_hold_all] = solve(psi_hold_all, mesh, matprop, pdecoef, bc); 
params.energy0 = 0.5 * dot(F,U); % initial energy 
max_vol = vol_hold_all(1); %volinit = volinit.*gamma; %calculate maximum vol for each fase
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
%figure(2); clf;
%set(2,'WindowStyle','docked');
pdeplot(p,e,t,'xydata',ones(size(p,2),1),'xystyle','flat','colormap','gray',...
              'xygrid','off','colorbar','off'); axis image; axis off; 
k = 1; iter = 0;
gsf = []; gEpot = []; gth = []; gvol = [];
option = 'null';

%while not(strcmp(option,'s'))
    %-----------Initial values-----------%
    [U,F,volume] = solve(psi, mesh, matprop, pdecoef, bc);
    [sf,Epot] = shfunc(F,U,volume,params);
    
    % compute function g: g = -DT (bulk) and g = DT (inclusion)
    dt = topder(mesh,U,volume,matprop,psi,params,pdecoef);
    dt = dt/normL2(unitM,dt);  % g function normalization
    
    dt(phold,:) = psi(phold,:); % freeze dt function
    
    for i = 1 : length(matprop.E)-1
        cosin(i) = psi(:,i)'*unitM*dt(:,i);
    end
    cosin = max(min(sum(cosin),1.0),-1.0);
    theta = max(real(acos(cosin)),1.0e-4);
    
    difvol = volume(1:end-1)-voltarget;
    ic = (abs(difvol) > volstop); %index control
    while and(and( or(any(ic),theta > stop) , k/2 > kmin), iter<=4) % remove iter<=4
                
        iter = iter + 1;
        [U,F,volume] = solve(psi, mesh, matprop, pdecoef, bc);
        [sf,Epot] = shfunc(F,U,volume,params);
        
        % compute function g: g = -DT (bulk) and g = DT (inclusion)
        dt = topder(mesh,U,volume,matprop,psi,params,pdecoef);
        dt = dt./normL2(unitM,dt);  % g function normalization
        
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
            [sf,Epot] = shfunc(F,U,volume,params);
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
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        
        facecolor = [186 212 244]/255;
        linecolor = [0 62 102]/255;
        %HH = figure(3); clf; %HH.Units='Normalized'; HH.OuterPosition=[0 0 1 1];
%        set(3,'WindowStyle','docked');
        multimat_plot( p,t,fi );
        drawnow
%             FRAMES(iter) = getframe(gcf);
%         figure(4); clf;hold on, plot(gsf,'k-','LineWidth',1.2); %title('Shape Function');
%         yyaxis left; ylabel('Shape Function'),
%         yyaxis right; plot(gEpot,'-','LineWidth',1.2);   ylabel('Epot'),grid minor, xlim([0 iter])
%         set(4,'WindowStyle','docked');
%         figure(5); clf; plot(gth*180/pi,'-','MarkerFaceColor',facecolor,'Color',linecolor,'LineWidth',1.2); title('Theta Angle');
%         grid minor , xlim([0 iter]),ylim([min(gth*180/pi)-5 max(gth*180/pi)+10])
%         set(5,'WindowStyle','docked');   
%         drawnow
%         figure(6); clf;
%         for index = 1:size(gvol,1)
%             plot(gvol(index,:),'-','LineWidth',1.2), hold on
%             leyenda{index} = sprintf('$V_{f%d}$',index);
%         end
%         title('Volume'); legend(leyenda),grid minor, xlim([0 iter]), ylim([-5 105])    
%         set(6,'WindowStyle','docked');
        
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
                penalty = penalty + (tau(ic)./auglag(ic)) .* (max(0,penalty - auglag(ic).*(1.0-coef))-penalty);
                params.penalty(ic) = penalty;
                k = 1;
            end
        end
    end
        
%     option = 'null';

%     if (option == 'r')
%         for i = 1:(strcmp(mesh.remesh,'longest')*2 + strcmp(mesh.remesh,'regular'))
%             [p,e,t,psi] = refinemesh(g,p,e,t,psi,mesh.remesh);
%             mesh.p = p; mesh.e = e; mesh.t = t;
%         end
%         
%         % update mesh parameters
% %         figure(1); clf;
% %         set(1,'WindowStyle','docked');
% %         pdemesh(p,e,t); axis image; axis off;
%         [~,unitM,~] = assema(p,t,0,1,0);
%         area = pdetrg(p,t); mesh.area = area;
%         k = 1;
%         
%         %Update bc
%         ldir = bc.Neu.ldir;
%         ldof = bc.Neu.ldof;
%         lval = bc.Neu.lval;
%         pbc = bouncon(p,e,t,ldir,ldof,lval);
%         pbc = unique(pbc,'rows'); %just for symetry_cog_2 eliminate otherwise
%         bc.pNeu = pbc;
%         
%         ldir = bc.Dir.ldir;
%         ldof = bc.Dir.ldof;
%         lval = bc.Dir.lval;
%         pbc = bouncon(p,e,t,ldir,ldof,lval);
%         pbc = unique(pbc,'rows'); %just for symetry_cog_2 eliminate otherwise
%         bc.pDir = pbc;
%         
%         % freeze groups
%         phold = [];
%         for i = 1:size(mesh.ghold,1)
%             phold = cat(1,phold,unique(t(1:3,t(4,:)==mesh.ghold(i)))); %nodes to block
%         end
%         pfree = setdiff(1:size(p,2),phold);
%         option = 'null';
%         
%         % shape functional and volume associated to the hold-all domain
%         psi_hold_all = ones(length(p),length(matprop.E)-1); psi_hold_all(:,1) = -1;
%         [U,F,vol_hold_all] = solve(psi_hold_all, mesh, matprop, pdecoef, bc);
%         params.energy0 = 0.5 * dot(F,U); % initial energy
%         max_vol = vol_hold_all(1);
%         params.volinit = max_vol; % store maximum volume for each fase
%         voltarget = max_vol.*volfrac;
%         
%         % stop criterion over the volume constraint
%         volstop = voleps.*max_vol;
%         if penalization == 1
%             volstop = max_vol;
%         end
%     end    
%end

end