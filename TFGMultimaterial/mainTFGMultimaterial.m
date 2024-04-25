function psi = mainTFGMultimaterial()


clear all; close all; clc
% load problem data

% Topology Optimization parameters: set penalization, line search and stop criterion parameters
params = ParametersComputer();
params = params.params;

% Generate mesh
mesh = MeshComputer();

% Generate boundary conditions
cParams.e = mesh.e;
cParams.p = mesh.p;
cParams.t = mesh.t;
cParams.g = mesh.g;
bc = BoundaryConditionsComputer(cParams);

% Set material properties
matProp = MaterialPropertiesComputer();

nMat = 4;
% Set pde coefficients
mat.A = matProp.matA;
mat.B = matProp.matB;
mat.C = matProp.matC;
mat.D = matProp.matD;
pdeCoeff = PDECoefficientsComputer(mat);


% Set an initial guess on level set
ls = InitialLevelSetComputer(mat,cParams);
psi = ls.psi;

% Plot mesh
figure(1); clf;
pdeplot(mesh.p,mesh.e,mesh.t,'xydata',mesh.t(4,:),'xystyle','flat','colormap','gray',...
              'xygrid','off','colorbar','off','title','Geometry'); 
axis image; axis off;

pdemesh(mesh.p,mesh.e,mesh.t); axis image; axis off;

pdegplot(mesh.g,"EdgeLabels","on","FaceLabels","on") % es per veure la geometria de la mesh

%freeze nodes --- WHY??
phold = [];
for i = 1:size(mesh.ghold,1)
    phold = cat(1,phold,unique(mesh.t(1:3,mesh.t(4,:)==mesh.ghold(i))));
end
pfree = setdiff(1:size(mesh.p,2),phold); % free nodes


s.connec = mesh.t';
s.connec = s.connec(:,1:3);
s.coord  = mesh.p';
%m = Mesh.create(s);
m = TriangleMesh(2,1,80,40); % why connectivites and coord are not the same? 


% shape functional and volume associated to the hold-all domain
psi_hold_all = ones(length(mesh.p),nMat-1);
psi_hold_all(:,1) = -1;
s.matProp = matProp;
s.mesh = mesh;
s.pdeCoeff = pdeCoeff;
s.psi = psi_hold_all;
s.bc  = bc.bc;


sys = FEMSolver(s);
[U,F] = sys.computeStiffnessMatrixAndForce();
s.tfi = sys.charFunc;
s.mesh = mesh;
vol = VolumeComputer(s);
vol_hold_all = vol.computeVolume();

% ElasticProblem
% [U,F]


params.energy0 = 0.5 * dot(F,U); % initial energy 
max_vol = vol_hold_all(1); %volinit = volinit.*gamma; %calculate maximum vol for each fase
params.max_vol = max_vol; % store maximum volume for each fase
voltarget = max_vol.*params.volfrac;

% stop criterion over the volume constraint 
volstop = params.voleps.*max_vol;
if params.penalization == 1
    volstop = max_vol;
end

p = mesh.p;
t = mesh.t;
e = mesh.e;



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
    %[U,F,volume] = solve(psi, mesh, matprop, pdecoef, bc);
        
    sys = FEMSolver(s);
    [U,F] = sys.computeStiffnessMatrixAndForce();
    s.tfi = sys.charFunc;
    s.mesh = mesh;
    vol = VolumeComputer(s);
    volume = vol.computeVolume();

    [sf,Epot] = shfunc(F,U,volume,params);
    
    % compute function g: g = -DT (bulk) and g = DT (inclusion)
    dt = topder(mesh,U,volume,matProp,psi,params,pdeCoeff);
    dt = dt/normL2(unitM,dt);  % g function normalization
    
    dt(phold,:) = psi(phold,:); % freeze dt function
    
    for i = 1 : (nMat-1)
        cosin(i) = psi(:,i)'*unitM*dt(:,i);
    end
    cosin = max(min(sum(cosin),1.0),-1.0);
    theta = max(real(acos(cosin)),1.0e-4);
    
    difvol = volume(1:end-1)-voltarget;
    ic = (abs(difvol) > volstop); %index control
    while and(and( or(any(ic),theta > params.stop) , k/2 > params.kmin), iter<=4) % remove iter<=4
    % while and( or(any(ic),theta > params.stop) , k/2 > params.kmin)           
        iter = iter + 1;

        s.psi = psi;
        s.mesh = mesh;
        s.matProp = matProp;
        s.pdeCoeff = pdeCoeff;
        s.bc  = bc.bc;

        sys = FEMSolver(s);
        [U,F] = sys.computeStiffnessMatrixAndForce();
        s.tfi = sys.charFunc;
        s.mesh = mesh;
        vol = VolumeComputer(s);
        volume = vol.computeVolume();
        [sf,Epot] = shfunc(F,U,volume,params);
        
        % compute function g: g = -DT (bulk) and g = DT (inclusion)
        dt = topder(mesh,U,volume,matProp,psi,params,pdeCoeff);
        dt = dt./normL2(unitM,dt);  % g function normalization
        
        dt(phold,:) = psi(phold,:); % freeze dt function
        
        for i = 1 : (nMat-1)
            cosin(i) = psi(:,i)'*unitM*dt(:,i);
        end
        cosin = max(min(sum(cosin),1.0),-1.0);
        theta = max(real(acos(cosin)),1.0e-4);
        
        % performs a line-search
        sfold = sf; psiold = psi; sf = sf + 1; k = min(1-(1e-5),1.5*k); % trick
%         
        while and((sf > sfold) , k>params.kmin)
%             (sf - sfold > eps(sf))
            % update level-set function
            psi(pfree,:)  = (sin((1-k)*theta)*psiold(pfree,:)...
                +  sin(k*theta)*dt(pfree,:))./sin(theta);
            
            s.psi = psi;
            s.mesh = mesh;
            s.matProp = matProp;
            s.pdeCoeff = pdeCoeff;
            s.bc  = bc.bc;

            sys = FEMSolver(s);
            [U,F] = sys.computeStiffnessMatrixAndForce();
            s.tfi = sys.charFunc;
            s.mesh = mesh;
            vol = VolumeComputer(s);
            volume = vol.computeVolume();
            [sf,Epot] = shfunc(F,U,volume,params);
            k = k / 2;
        end
        psi = psi/normL2(unitM,psi);
        k = k * 2;

        s.psi = psi;
        s.p = mesh.p;
        s.t = mesh.t;

        charfun = CharacteristicFunctionComputer(s); % s'ha de construir la classe - charfunc!!
        [fi,~] = charfun.compute();
        %[fi,~] = charfunc( p,t,psi); %calculate characteristic function for plot
        
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
        %     FRAMES(iter) = getframe(gcf);
        % figure(4); clf;hold on, plot(gsf,'k-','LineWidth',1.2); %title('Shape Function');
        % yyaxis left; ylabel('Shape Function'),
        % yyaxis right; plot(gEpot,'-','LineWidth',1.2);   ylabel('Epot'),grid minor, xlim([0 iter])
        % set(4,'WindowStyle','docked');
        % figure(5); clf; plot(gth*180/pi,'-','MarkerFaceColor',facecolor,'Color',linecolor,'LineWidth',1.2); title('Theta Angle');
        % grid minor , xlim([0 iter]),ylim([min(gth*180/pi)-5 max(gth*180/pi)+10])
        % set(5,'WindowStyle','docked');   
        % drawnow
        % figure(6); clf;
        % for index = 1:size(gvol,1)
        %     plot(gvol(index,:),'-','LineWidth',1.2), hold on
        %     leyenda{index} = sprintf('$V_{f%d}$',index);
        % end
        % title('Volume'); legend(leyenda),grid minor, xlim([0 iter]), ylim([-5 105])    
        % set(6,'WindowStyle','docked');
        
        %Update augmented lagrangian parameters
        difvol = volume(1:end-1)-voltarget;
        ic = (abs(difvol) > volstop); %index control
        
        if any(ic==1) 
            if params.penalization == 2 % increase the penalization parameter
                params.penalty = 2.0 * params.penalty;
                k = 1;
            elseif params.penalization == 3 % increase the lagrangian multiplier
                coef = volume(1:end-1)./ voltarget; coef = coef(ic); tau = params.auglag;
                penalty = params.penalty(ic);
                penalty = penalty + (tau(ic)./params.auglag(ic)) .* (max(0,penalty - params.auglag(ic).*(1.0-coef))-penalty);
                params.penalty(ic) = penalty;
                k = 1;
            end
        end
    end
        
end