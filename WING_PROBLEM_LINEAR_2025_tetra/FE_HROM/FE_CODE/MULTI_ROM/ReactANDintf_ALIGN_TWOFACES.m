function [DATAOUT,BASES,BasisINT,TEXTP,BasisUdef] =  ReactANDintf_ALIGN_TWOFACES(DATA_REFMESH,...
    DATAOUT,DATAIN,BasisUdef,BasisRdef,...
    SinvVal_Udef,f,fI,SinvVal_Rdef,FACES_GROUPS,BASES,TEXTP)

if nargin == 0
    load('tmp2.mat')
end

%error('Unreliable method-......')
% Adaptation of ALIGNMENT_METHOD.m

% ROTATION MATRICES
DATAIN.RotationMatrixLocal = DATA_REFMESH.RotationMatrixFace ;

ndim = size(DATA_REFMESH.COOR,2);

% Rigid body modes for interface
Vrb = DATA_REFMESH.RigidBodyModesInterface ;
DATAOUT.BasisIntRB = Vrb ;

% STEP 1) DISPLACEMENT MODES (ALIGNED WITH THE REACTION MODES )
% -------------------------------------------------------------
 
[BasisUdef,SinvVal_Udef, TEXTP]= ReactionDispModesEqual(BasisUdef,BasisRdef,f,...
    SinvVal_Udef,TEXTP) ;
DATAOUT.BasisUdef = BasisUdef ;
BASES.DISPLACEMENTS.U = BasisUdef ;
BASES.DISPLACEMENTS.S = SinvVal_Udef ;
DATAOUT.BasisUdef = BasisUdef; 
DATAOUT.SingVal_Udef = SinvVal_Udef; 
% PLOT AGAING DISPL. MODES
BasisUdefPlot = BasisUdef ;
CNref = DATA_REFMESH.CN  ;
COOR = DATA_REFMESH.COOR  ;
NAME_MODES = [DATAIN.NAME_WS_MODES(1:end-4),'DISP_ALIGNED' ];
DATA= [] ;
LEGENDG = 'DISP.'
TEXTP = GidPostProcessModes_dom(COOR,CNref,DATA_REFMESH.TypeElement,BasisUdefPlot,DATA_REFMESH.posgp,...
    NAME_MODES,DATA,LEGENDG,TEXTP);
DATAOUT.BasisRdef = BasisRdef ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------
% Geometric mass matrices
% ------------------------
nfaces = length(Vrb);
M1d = DATA_REFMESH.GeometricMassMatrixInterface ;
M = cell(nfaces,1) ;
for iface= 1:nfaces
    M{iface} = sparse(size(M1d{iface},1)*ndim,size(M1d{iface},2)*ndim) ;
    for idim = 1:ndim
        M{iface}(idim:ndim:end,idim:ndim:end) = M1d{iface} ;
    end
end
DATAOUT.MassMinterfaces = M ;


% -----------------------------------------------------------------------
% DETERMINATION OF INTERFACE MODES.
% ----------------------------------
Mdom = DATA_REFMESH.M ;
BasisUrb =   DATA_REFMESH.BasisUrb ;
% Candidates FACE,BY,FACE  (AFTER COMPUTING THE INTERSECTION)
% ----------------------
[BasisINTFdef,BasisINTrbORTH,BasisRdefROT,TEXTP,BasisUdomROT,BasisINTFdef_VIRGIN  ]= ...
    CandidatesInterfaceModesSTD(DATAIN,Vrb,DATAOUT,...
    nfaces,fI,BasisUdef,M,...
    TEXTP,FACES_GROUPS,SinvVal_Udef,BasisUrb,Mdom,BasisRdef,TEXTP,SinvVal_Rdef) ;

%% Rotation ---- Intersection   (Commented out 20-Jan-2020)
%BasisRdef1 = [BasisRdefROT{1}; BasisRdefROT{2}] ;
%BasisRdef2 = [BasisRdefROT{3}; BasisRdefROT{4}] ;


%%%%% % Commented out.... 20-Jan.2020

% SinvVal_Rdef  = SinvVal_Rdef./SinvVal_Rdef(1) ;
% 
% DATAIN= DefaultField(DATAIN,'INCLUDE_SINGULAR_VALUES_COMPUTATION_OF_INTERSECTIONS',0) ; %
% if DATAIN.INCLUDE_SINGULAR_VALUES_COMPUTATION_OF_INTERSECTIONS ==1
%     BasisRdef1  = bsxfun(@times,BasisRdef1',SinvVal_Rdef(1:size(BasisRdef1,2)) )' ;
%     BasisRdef2   = bsxfun(@times,BasisRdef2',SinvVal_Rdef(1:size(BasisRdef1,2)))' ;
% end
%%%%%


%%%%% % Commented out.... 20-Jan.2020

% %% Orthogonal basis matrix for the intersection between BasisRdef1 and BasisRdef2  ---> R
% DATAIN= DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_AlignmentMethod',1e-4) ; %
% DATAIN= DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_ANGLES_INTERSECTION_RU',1) ; % =
% DATAIN= DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_ANGLES_INTERSECTION_V',0.001) ; % =
% 
% TOL_ANGLE =DATAIN.TOL_DeformationalInterfaceModes_ANGLES_INTERSECTION_RU;
% TOL = DATAIN.TOL_DeformationalInterfaceModes_AlignmentMethod;
% TEXTP{end+1} = '------------------------------------------------------------';
% TEXTP{end+1} = ['2-FACES reactions'];
% TEXTP{end+1} = '------------------------------------------------------------';
% 
% % I think this is redundant (3-Oct-2019).... 14-Jan-2020... Why  ? 
% [R,ANGLES,sA,sB] = IntersectionSubspaces(BasisRdef1,BasisRdef2,TOL,TOL_ANGLE);
% 
% TEXTP{end+1} = ['Number of reaction modes  = ', num2str(size(BasisRdef1,2)) ];
% TEXTP{end+1} = ['Number of reaction modes after truncation FACE -  = ', num2str(sA) ];
% TEXTP{end+1} = ['Number of reaction modes after truncation FACE +  = ', num2str(sB) ];
% TEXTP{end+1} = ['Dimension of reactions intersection space = ', num2str(size(R,2)) ];
% TEXTP{end+1}  = ['---------------------------------------------'] ;


%nmodesR_12 = size(R,2) ;
nRB  =size(Vrb{1},2) ;
nmodesTOTAL = size(BasisUdef,2) + nRB; 

nmodesU1 = nRB + size(BasisINTFdef{1},2) ;
nmodesU2 = nRB + size(BasisINTFdef{2},2) ;
% Check wether this criterion makes sense.... 
if 2*(nmodesU1 + nmodesU2) <= nmodesTOTAL  % 2*(nmodesU1 + nmodesU2)  <= size(BasisRdef,2) + nRB 
    % 18th-Sept-2019: What's the usefulness of this criterion ? 
 
    
    TEXTP{end+1} =['No need to use kinematic constraints.....'] ;
    BasisINT{1} = [Vrb{1},BasisINTFdef{1}] ;
    BasisINT{2} = [Vrb{2},BasisINTFdef{2}] ;
    BasisINT{3} = BasisINT{1}   ;
    BasisINT{4} = BasisINT{2}   ;
    %
    %
    DATAOUT.BasisInt = BasisINT ;
    %
    
else
    % error('Not imp. ')
    % KINEMATIC APPROACH IS REQUIRED
    
    DATAIN = DefaultField(DATAIN,'AlignmentMethod_withKINEMATIC_USE_VIRGIN_MODES',1) ;     
    if DATAIN.AlignmentMethod_withKINEMATIC_USE_VIRGIN_MODES == 1
      BasisINTFdef =   BasisINTFdef_VIRGIN ; 
    else
        error('This option proved to be unreliable  ')
      %  DmatrixImplementation.tex
    end
    
     [TEXTP,DATAOUT,BasisINT]  =KinemApproachCandidatsDOFS_LOC(BasisINTFdef,Vrb,fI,BasisUrb,BasisUdef,...
    Mdom,TEXTP,DATAOUT,BasisINTrbORTH,BasisUdomROT,DATAIN,M) ; 
    
end