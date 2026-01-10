function [BasisRdef,BasisINT,MSG,BasisUdef,SingVal_Udef,Vdef_Morth] = ...
    ReactionAndInterfaceLocalModes_UNIFIED(BasisUdef,BasisRdef,f1,f2,TOL_SINGULAR_VALUES_Hqr,...
    nBASES_BEAM,DATA_REFMESH,Vrb,M,DATAIN,SingVal_Udef,SingVal_Rdef,BasisUrb,Mdom,MSG)
if nargin == 0
    load('tmp2.mat')
end

f = [f1;f2] ;

% Candidates to be interface displacement modes
% ----------------------------------------------
% Rotation
% ---------
if ~isempty(DATA_REFMESH.RotationMatrixFace)
    iface = 2;
    DATAIN.ROTATION_LOC = DATA_REFMESH.RotationMatrixFace{iface} ;
    R = DATAIN.ROTATION_LOC ;
else
    R = [] ;
end


% Selection of reaction modes
% ---------------------------
DATAIN = DefaultField(DATAIN,'EQUAL_NUMBER_REACTIONS_NUMBER_MODES_REACTION_DRIVEN',1) ; %Method, 26-Feb-2019


DATAIN = DefaultField(DATAIN,'TOL_SINGULAR_VALUES_Hqr_REACTIONS',1e-4) ;
nrigidB = size(Vrb,2) ;

if  DATAIN.EQUAL_NUMBER_REACTIONS_NUMBER_MODES_REACTION_DRIVEN == 0
    error('Method abandoned in August 2019')
    [BasisRdef, TEXTR]= ReactionModesLocalEffectsNEW(BasisUdef,BasisRdef,f,...
        DATAIN.TOL_SINGULAR_VALUES_Hqr_REACTIONS,SingVal_Rdef,nrigidB,DATAIN) ;
else
    DATAIN = DefaultField(DATAIN,'KINEMATIC_CONSTRAINTS_MODES',[]);   % NEW METHOD, 14-MAY-2019
    
    DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'DEIM_BASED_METHOD',0);   % NEW METHOD, 14-MAY-2019
    
    if DATAIN.KINEMATIC_CONSTRAINTS_MODES.DEIM_BASED_METHOD == 1
        %     if ndim == 2
        %         nRB = 3;  % Number of rigid body modes
        %     else
        %         nRB  = 6 ;
        %     end
        nTOTAL =  nrigidB + DATAIN.nMODES_REACTIONS ;
        if mod(nTOTAL,2) ~=0
            DATAIN.nMODES_REACTIONS = DATAIN.nMODES_REACTIONS-1;
        end
        BasisRdef = BasisRdef(:,1:DATAIN.nMODES_REACTIONS)  ;
    end
    
    
    
    % Determination of aligned displacement method
    % Modified  on May-13th-2020, to filter out modes that  do not
    % contribute to the work done by the  unit cell
    [BasisUdef,SingVal_Udef, MSG,BasisRdef,SingularValuesNewREAC]= ...
        ReactionDispModesEqual(BasisUdef,BasisRdef,f,...
        SingVal_Udef,MSG,DATAIN) ;
    
    if ~isempty(SingularValuesNewREAC)
        SingVal_Rdef = SingularValuesNewREAC ; 
    end
    
    
    
    % PLOT AGAING DISPL. MODES
    
    BasisUdefPlot = BasisUdef ;
    
    % end
    
    CNref = DATA_REFMESH.CN  ;
    COOR = DATA_REFMESH.COOR  ;
    NAME_MODES = [DATAIN.NAME_WS_MODES(1:end-4),'DISP_ALIGNED' ];
    DATA= [] ;
    LEGENDG = 'DISP.' ;
    
    MSG{end+1} = '-----------------------------------';
    MSG{end+1} = 'DISPLACEMENT MODES AFTER ALIGNMENT';
    MSG{end+1} = '-----------------------------------';   
    MSG =   GidPostProcessModes_dom(COOR,CNref,DATA_REFMESH.TypeElement,BasisUdefPlot,DATA_REFMESH.posgp,...
        NAME_MODES,DATA,LEGENDG,MSG);   
    MSG{end+1} = '-----------------------------------';
    
    
    
    
end





%DATAOUT.BasisRdef = BasisRdef ;


% Candidate for being interface modes
% Vdef_Morth --> Rigid body modes (always included)
% Vrv_Morth --> Additional modes
[Vdef_Morth,Vrb_Morth] = CandidatesInterfaceModes_UNIFIED(BasisUdef,f1,f2,SingVal_Udef,...
    DATAIN,M,Vrb,BasisRdef,BasisUrb) ;

% warning('borrar')
% Vcand = [Vrb_Morth,Vdef_Morth] ;
% T_1 = BasisRdef(f1,:)'*Vcand ;
% T_2 = BasisRdef(f2,:)'*Vcand ;



DATAIN = DefaultField(DATAIN,'DeformationalInterfaceModes_AlignmentMethod',1) ;
DATAIN = DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_AlignmentMethod',1e-3) ;


%DATAIN.KINEMATIC_CONSTRAINTS_MODES.DEIM_BASED_METHOD = 1;   % NEW METHOD, 14-MAY-2019
DATAIN = DefaultField(DATAIN,'KINEMATIC_CONSTRAINTS_MODES',[]) ;
DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'DEIM_BASED_METHOD',1) ;


if DATAIN.KINEMATIC_CONSTRAINTS_MODES.DEIM_BASED_METHOD == 1
    % May-14th- New method. Similar to the one implemented for RVEs. See
    % InterfaceRVEKinematicConstraint.m
    error('Method abandoned in Aug-2019 in favor of DATAIN.KINEMATIC_CONSTRAINTS_MODES.DEIM_BASED_METHOD ==0 ')
    BasisUdom = [BasisUrb,BasisUdef] ;
    
    [BasisINT,TEXTR ]= InterfaceRVEKinematicConstraint_BEAMS(Vrb,M,BasisRdef,Vdef_Morth,...
        TOL_SINGULAR_VALUES_Hqr,f1,f2,R,DATAIN,Vrb_Morth,Mdom,BasisUdom,BasisUrb,TEXTR) ;
    
    
else
    
    if DATAIN.DeformationalInterfaceModes_AlignmentMethod == 1
        % warning('This implementation is not reliable')
        DATAIN.RotationMatrixLocal{1} = [] ;
        DATAIN.RotationMatrixLocal{2} = R' ;
        iface = 1;
        [BasisINT,~,MSG ]=  DeformModesInterface_AlignmentMethod(BasisRdef,f1,f2,...
            Vrb,M,DATAIN,BasisUdef,SingVal_Udef,SingVal_Rdef,iface,MSG,Vdef_Morth) ;
        % TEXTB ={} ;
    else
        error('Method abandoned in Aug-2019 ')
        % --------
        %% Final selection of interface modes
        [BasisINT,TEXTR ]= SelectCandidatesInterfaceModes_BEAMSnew(Vrb,M,BasisRdef,Vdef_Morth,...
            TOL_SINGULAR_VALUES_Hqr,f1,f2,R,DATAIN,Vrb_Morth) ;
        
    end
    
end
