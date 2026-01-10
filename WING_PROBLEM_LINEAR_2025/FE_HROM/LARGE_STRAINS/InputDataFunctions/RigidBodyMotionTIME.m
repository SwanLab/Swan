function RIGID_BODY_MOTION = RigidBodyMotionTIME(DATA,ndim,MESH,GEOproperties,DATALOC)
% Computation of uRB  for all time steps
% See 07_RotatingFrameStatic.pdf
% JAHO 7-Jan-2021 (Borovet,Bulgaria)
% ----------------------------------------

if nargin == 0
    load('tmp1.mat')
end


% Global rotation  (See typical input data file
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/
% HROM/07_BEAM_rotation_specialFORM/DirichletLocal.m)

% ---------------
% iloadstate_loc = iloadstate;
% ROTATION_GLO = [] ;
% ROTATION_GLO(iloadstate_loc).CENTER.COOR = [] ;
% ROTATION_GLO(iloadstate_loc).CENTER.ISCENTROID = 1 ;
% ROTATION_GLO(iloadstate_loc).Z.MAXANGLE = rotation_GLO;
% ROTATION_GLO(iloadstate_loc).Z.TIMEFUN = @(t) t/tEND ;
% ROTATION_GLO(iloadstate_loc).INTERVAL  = [0,tEND] ;
%
% TRANSLATION_GLO = [] ;
% TRANSLATION_GLO(iloadstate_loc).TIMEFUN = @(t) t/tEND ;
% TRANSLATION_GLO(iloadstate_loc).AMPLITUDE = [disp_x_face1,disp_y_face1] ;
% TRANSLATION_GLO(iloadstate_loc).INTERVAL  = [0,tEND] ;
%
% RB_MOTION.TRANSLATION  = TRANSLATION_GLO ;
% RB_MOTION.ROTATION  = ROTATION_GLO ;


%
TRANSLATION = DATALOC.RB_MOTION.TRANSLATION ;
ROTATION= DATALOC.RB_MOTION.ROTATION ;
%
COOR = MESH.COOR  ;
nloads = length(TRANSLATION);
ntimesteps = length(DATA.STEPS) ;
U = cell(1,nloads) ;
a = cell(nloads,1) ;
RotationMatrix = cell(1,ntimesteps) ; 

for  iload = 1:nloads
    INTERVAL =  TRANSLATION(iload).INTERVAL ;
    
    ROTATION_GLO = ROTATION(iload) ;
    TRANSLATION_GLO = TRANSLATION(iload) ;
    
    % Rotation center for this loading stage
    if isempty(ROTATION_GLO.CENTER.COOR)
        isurfGLO =  ROTATION_GLO.CENTER.ISCENTROID  ;
        CENTERglo  =  GEOproperties.FACES{isurfGLO}.CENTROID;
    else
        CENTERglo =ROTATION_GLO.CENTER.COOR ;
    end
    ndof = prod(size(COOR));
    D= zeros(ndof,ntimesteps)   ;
    Dloc = zeros(ndof,1) ; dLOCAL = zeros(ndim,size(COOR,1)) ;
    for istep = 1:ntimesteps
        
        TIMELOC = DATA.STEPS(istep) ;
        
        %
        %         Dloc =  D(:,istep)  ;
        %         % RIGID BODY MOTION RELATIVE TO THE CENTROID OF THE FACE UNDER
        %         % CONSIDERATION
        %         Dloc = RigidBodyMotionALL(TRANSLATION,INTERVAL,TIMELOC,rnodLOC,Dloc,ndim,ROTATION,COORfaceREL) ;
        % D(:,istep)  = Dloc ;
        %
        
        %% GLOBAL ROTATION
        % --------------------
        %   if ~isempty(ROTATION_GLO)
        % Global rotation
        % ------------------
        %     dLOCAL = D(:,istep) ; %
        %    Dloc = reshape(dLOCAL,ndim,[]) ;
        %   COORface_NEW = dLOCAL + COORface' ; %
        COORrel = bsxfun(@minus,COOR',CENTERglo')' ;
        %  Dloc = zeros(size(Dloc));
        %  else
        %     COORfaceREL_new = [] ;
        % end
        %    if ~isempty(ROTATION_GLO) ||      ~isempty(TRANSLATION_GLO)
        rnodLOC = 1:size(COOR,1) ;
        [D(:,istep),RotationMatrix_step]= RigidBodyMotionALL_compound(TRANSLATION_GLO,INTERVAL,TIMELOC,rnodLOC,Dloc,ndim,ROTATION_GLO,COORrel,dLOCAL) ;
        
        %   end
        
        if any(RotationMatrix_step)
            RotationMatrix{istep} = RotationMatrix_step ; 
        end
        
        
        
        
        
    end
    
    % Next we apply the SVD
    if any(any(D))
        
        % What if D does not fit into memory ... ? 
        % What if the SVD is too slow ? 
        SIZE_D = (8e-6)*prod(size(D)) ; 
        LIMIT = 50 ; % Mbytes 
        
        nclusters = ceil(SIZE_D/LIMIT) ;
        
        ssss = [] ; 
       NSTEPS_size= [] ;
        % We have to divide  1:length(DATA.STEPS) into nclusters
        FREQ  = ceil(size(D,2)/nclusters) ;
        iini = 1;
        %ifin = FREQ ;
        for i=1:nclusters
            ifin = iini + FREQ-1 ;
            ifin = min(ifin,size(D,2)) ;
            NSTEPS_size =[NSTEPS_size,  (ifin-iini+1)] ;
            iini = ifin +1 ;
        end
        
        
        if length(NSTEPS_size) == 1
            [Uloc,SS,VV] =  RSVDT(D) ;
        else
            D = mat2cell(D,size(D,1),NSTEPS_size) ; 
            epsilonD = 0*ones(length(D)) ; 
            [Uloc,SS,VV] =  RSVDqp(D,epsilonD) ;
        end
        
        
        
        aloc = bsxfun(@times,VV',SS) ;
        U{iload} =  Uloc ;
        a{iload} =  aloc ;
    else
        
        U{iload}  = zeros(size(D,1),1) ;
        a{iload} = zeros(1,size(D,2)) ;
    end
    
    
    %     U{izone}(:,iload) = displ{iload} ;    % Uniform Boundary Conditions
    %     FactorSteps =  LOCVAR.PRESCRIBED_DISP(iload).TIMEFUN(DATA.STEPS) ;
    %     LIMITS_interval =  LOCVAR.PRESCRIBED_DISP(iload).INTERVAL ;
    %     FactorSteps = FactorSteps.*(DATA.STEPS >= LIMITS_interval(1)) ;
    %     a{izone}(iload,:) = FactorSteps.*(DATA.STEPS <= LIMITS_interval(2)) ;
end


U = cell2mat(U) ;
a = cell2mat(a) ;


RIGID_BODY_MOTION.U = U ; 
RIGID_BODY_MOTION.a = a ; 
RIGID_BODY_MOTION.RotREFP = RotationMatrix ; 




