function [DOFr,dR] = DirichletCONDtimeCABLE(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC)
% Goal. Determine DOFr and    dR(t)   (prescribed DOFs and value of the corresponding
% displacements)  in a space-time separated fashion
%
%
% JAHO- 23-JUN-2O22
%--------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
elseif nargin == 5
    DATALOC = [] ;
end

nzones = length(DIRICHLET) ; % number of zones with prescribed DOFs
DOFr = cell(nzones,1); U = cell(nzones,1) ; a = cell(nzones,1) ;
for  izone = 1:length(DIRICHLET)
    
    LOCVAR=  DIRICHLET(izone) ;
    
    %%%%
         rnodLOC=  MESH.POINTS_BOUNDARY(DIRICHLET(izone).NUMBER_POINT) ; 
    
    DOFloc = small2large(rnodLOC,ndim) ;  % Candidates for being DOFr
    
    % Determining the DOFr corresponding to this zone
    % Notice that the DOFr must be the same for all load states
    % ----------------------------------------------------------
    DOFrLOC =  cell(1,length(DIRICHLET(izone).PRESCRIBED_DISP)) ;
    displ =  cell(1,length(DIRICHLET(izone).PRESCRIBED_DISP)) ;
    nloads=  length(DIRICHLET(izone).PRESCRIBED_DISP) ;
    for iload = 1:nloads
        PRESCRIBED_DISPLACEMENT = DIRICHLET(izone).PRESCRIBED_DISP(iload);
        for idim = 1:ndim
            if ~isempty(PRESCRIBED_DISPLACEMENT.AMPLITUDE{idim})
                DOFloc_dim = DOFloc(idim:ndim:end) ;
                dddd =[ PRESCRIBED_DISPLACEMENT.AMPLITUDE{idim} ]*ones(length(DOFloc_dim),1) ;
                displ{iload} = [ displ{iload}; dddd] ;
                DOFrLOC{iload} = [DOFrLOC{iload};DOFloc_dim] ;  % DOFr corresponding
            end
        end
        if iload >1
            % Checking that the DOFs are the same
            if ~isempty(setdiff(DOFrLOC{iload},DOFrLOC{iload-1}))
                error('Ill-defined Boundary Conditions... Load states with distinct DOFs')
            end
        end
    end
    
    % Checking that there are no repeated indexes  (and remove them if indeed there are repetitions)
    if izone >1
        [REPEATED,iOLD,iNEW] = intersect(cell2mat(DOFr(1:izone-1)),DOFrLOC{1}) ;
        if ~isempty(REPEATED)
            IND =  setdiff(1:length(DOFrLOC{1}),iNEW) ;
            DOFrLOC = DOFrLOC{1}(IND) ;
            for iload =1:nloads
                displ{iload} = displ{iload}(IND) ;
            end
        else
            DOFrLOC = DOFrLOC{1} ;
        end
    else
        DOFrLOC = DOFrLOC{1} ;
    end
    %
    DOFr{izone} = DOFrLOC;
    U{izone} = zeros(length(DOFrLOC),nloads) ;
    a{izone} = zeros(nloads,length(DATA.STEPS)) ;
    
    
    for  iload = 1:nloads
        U{izone}(:,iload) = displ{iload} ;    % Uniform Boundary Conditions
        FactorSteps =  DIRICHLET(izone).PRESCRIBED_DISP(iload).TIMEFUN(DATA.STEPS) ;
        LIMITS_interval =  DIRICHLET(izone).PRESCRIBED_DISP(iload).INTERVAL ;
        FactorSteps = FactorSteps.*(DATA.STEPS >= LIMITS_interval(1)) ;
        a{izone}(iload,:) = FactorSteps.*(DATA.STEPS <= LIMITS_interval(2)) ;
    end
    U{izone} = sparse( U{izone}) ;
    
    % end
    
    
end

%dR.U = U ;
%dR.a = a ;

DOFr  =cell2mat(DOFr) ;
[DOFr,III ]= sort(DOFr);

dR.U = blkdiag(U{:}) ;
dR.a = cell2mat(a) ;

dR.U = dR.U(III,:) ;


%%% RIGID BODY MOTION
% -----------------------------
% Compute, for each time step
% 1) Translation  (3xtime_steps)
% 2) Relative coordinates
% 3) Rotation (3x3xtime_steps)

% DATALOC = DefaultField(DATALOC,'RB_MOTION',[]) ;
% 
% if ~isempty(DATALOC.RB_MOTION)
%     RIGID_BODY_MOTION = RigidBodyMotionTIME(DATA,ndim,MESH,GEOproperties,DATALOC)   ;
%     
%     % We have to redefine dR.U and dR.a ...
%     D = zeros(size(dR.U,1),size(dR.a,2)) ;
%     for itimestep = 1:size(dR.a,2)
%         uBAR = dR.U*dR.a(:,itimestep) ; % Original prescribed displacement
%         uBAR_rb  = RIGID_BODY_MOTION.U(DOFr,:)*RIGID_BODY_MOTION.a(:,itimestep) ; % Rigid body motion
%         uBARnew = (uBAR-uBAR_rb) ; %
%         uBARnew = reshape(uBARnew,ndim,[]) ;
%         uBARnew = RIGID_BODY_MOTION.RotREFP{itimestep}'*uBARnew ;
%         D(:,itimestep) =uBARnew(:);
%         
%         %
%     end
%     disp('Recomputing dR.U and dR.a')
%     TOLrel = 1e-10;
%     [U,S,V] = RSVDTrel(D,TOLrel) ;
%     dR.U = U ;
%     a = bsxfun(@times,V',S) ;
%     dR.a = a;
%     
%     
% else
%     RIGID_BODY_MOTION = [] ;
% end
% 
% dR.RIGID_BODY_MOTION = RIGID_BODY_MOTION;

