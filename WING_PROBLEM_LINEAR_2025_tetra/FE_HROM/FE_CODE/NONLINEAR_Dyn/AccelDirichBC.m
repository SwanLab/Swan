function [d0,v0,gBOUNDdd,gBOUNDd]=AccelDirichBC(DynDATA,ndofTOT,DOFm,DOFs,DOFl,Gbound,gBOUND,VECTOR_TIME)


% Acceleration  of Dirichlet boundary conditions
% First we address the issue of initial displacements and

if nargin == 0
    load('tmp.mat')
end

% velocities
switch DynDATA.d0
    case 'zero'
        d0 = zeros(ndofTOT,1) ;
    otherwise
        error('Option not implemented')
end
switch DynDATA.v0
    case 'zero'
        v0 = zeros(ndofTOT,1) ;
    otherwise
        error('Option not implemented')
end
% Newmark integration scheme
betaNM = DynDATA.betaNM ;
gammaNM = DynDATA.gammaNM ;

if ~isempty(DOFm)
    %  gvN = v0(DOFs) -Gbound*v0(DOFm) ;
    gN = d0(DOFs) -Gbound*d0(DOFm) ;
else
    %  gvN = v0(DOFs) ;
    gN = d0(DOFs)  ;
end


NEWMARK =0 ;
if  NEWMARK == 0
    
    gBOUNDdd = zeros(size(gBOUND));
    gBOUNDd  = zeros(size(gBOUND));
    istep=2 ;
    dt = VECTOR_TIME(istep)-VECTOR_TIME(istep-1);
    gNP1 = gBOUND(:,istep) ;
    gN  = gBOUND(:,istep-1) ;
    vNP1   = (gNP1-gN)/dt ;
    vN = vNP1 ;
    aN = zeros(size(vN));
    gBOUNDd(:,1) = vN ;
    
    for istep = 2:length(VECTOR_TIME)
        dt = VECTOR_TIME(istep)-VECTOR_TIME(istep-1);
        gNP1 = gBOUND(:,istep);
        vNP1 = (gNP1-gN)/dt ;
        aNP1 = (vNP1-vN)/dt ;
        %
        gBOUNDd(:,istep) = vNP1 ;
        gBOUNDdd(:,istep) = aNP1 ;
        vN = vNP1 ;
        aN = aNP1 ;
        gN = gNP1 ;
    end
    
    
else
    gvN = (gBOUND(:,2)-gBOUND(:,1))/(VECTOR_TIME(2)-VECTOR_TIME(1)) ;
    
    gaN = zeros(size(gN));
    gBOUNDdd = zeros(size(gBOUND));
    gBOUNDd  = zeros(size(gBOUND));
    gBOUNDd(:,1)  = gvN ;
    for istep = 2:length(VECTOR_TIME)
        gNP1 = gBOUND(:,istep);
        dt = VECTOR_TIME(istep)-VECTOR_TIME(istep-1);
        % Tilde displacement
        gNP1nm =  gN + dt*gvN + 0.5*(dt^2)*(1-2*betaNM)*gaN ;
        % Acceleration
        gaNP1 = (gNP1-gNP1nm)/dt^2/betaNM;
        % Tilde velocity
        gvNP1nm = gvN + (1-gammaNM)*dt*gaN ;
        % Velocity
        gvNP1 = gvNP1nm + gammaNM*dt*gaNP1 ;
        
        %% Updates
        gvN = gvNP1 ;
        gaN = gaNP1 
        gN  = gNP1  ;
        gBOUNDdd(:,istep) = gaNP1 ;
        gBOUNDd(:,istep) = gvNP1 ;
    end
    
end
%
% 
% inode = 10 ;
% figure(809)
% hold
% xlabel('Time (s)')
% ylabel('d (m)')
% 
% plot(VECTOR_TIME,gBOUND(inode,:),'r')
% 
% plot(VECTOR_TIME,gBOUNDd(inode,:),'g')
% plot(VECTOR_TIME,gBOUNDdd(inode,:),'k')