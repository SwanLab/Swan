function [DOFr,dR] = DirichletCONDtime_rigid(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC)
% Goal. Determine DOFr and    dR(t)   (prescribed DOFs and value of the corresponding
% displacements)  in a space-time separated fashion
%
% JAHO- 2-aPRIL-2O25,   HGs, Pedralbes,   Barcelona
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
        [Uloc,aloc,DOFrLOCC] =    DirichletCONDtimeRIGIDbody(LOCVAR,DATA,ndim,MESH,GEOproperties) ;
        U{izone} = sparse(Uloc) ;
        a{izone} = aloc ;
        DOFr{izone} = DOFrLOCC ; 
    
end

%dR.U = U ;
%dR.a = a ;

DOFr  =cell2mat(DOFr) ;
[DOFr,III ]= sort(DOFr);

dR.U = blkdiag(U{:}) ;
dR.a = cell2mat(a) ;

dR.U = dR.U(III,:) ;

% 
% %%% RIGID BODY MOTION
% % -----------------------------
% % Compute, for each time step
% % 1) Translation  (3xtime_steps)
% % 2) Relative coordinates
% % 3) Rotation (3x3xtime_steps)
% 
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
% 
