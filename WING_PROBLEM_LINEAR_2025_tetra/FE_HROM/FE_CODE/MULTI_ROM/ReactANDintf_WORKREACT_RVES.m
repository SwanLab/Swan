function [DATAOUT,BASES,BasisINT,TEXTP,BasisUdef] =  ReactANDintf_WORKREACT_RVES(DATA_REFMESH,DATAOUT,DATAIN,BasisUdef,BasisRdef,...
    SinvVal_Udef,f,fI,SinvVal_Rdef,FACES_GROUPS,BASES)

if nargin == 0
    load('tmp1.mat')
end
% Adaptation of ReactANDintf_WORKREACT.m  

ndim = size(DATA_REFMESH.COOR,2);

% Rigid body modes for interface
Vrb = DATA_REFMESH.RigidBodyModesInterface ;
DATAOUT.BasisIntRB = Vrb ;

% STEP 1) DISPLACEMENT MODES (ALIGNED WITH THE REACTION MODES )
% -------------------------------------------------------------
[BasisUdef,SinvVal_Udef, TEXTP]= ReactionDispModesEqual(BasisUdef,BasisRdef,f,...
    SinvVal_Udef) ;
DATAOUT.BasisUdef = BasisUdef ;
BASES.DISPLACEMENTS.U = BasisUdef ;
BASES.DISPLACEMENTS.S = SinvVal_Udef ;

% PLOT AGAING DISPL. MODES
BasisUdefPlot = BasisUdef ;
CNref = DATA_REFMESH.CN  ;
COOR = DATA_REFMESH.COOR  ;
NAME_MODES = [DATAIN.NAME_WS_MODES(1:end-4),'DISP_ALIGNED' ];
DATA= [] ;
LEGENDG = 'DISP.'
GidPostProcessModes_dom(COOR,CNref,DATA_REFMESH.TypeElement,BasisUdefPlot,DATA_REFMESH.posgp,...
    NAME_MODES,DATA,LEGENDG);
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



 


% DETERMINATION OF INTERFACE MODES.
% ----------------------------------
Mdom = DATA_REFMESH.M ;
BasisUrb =   DATA_REFMESH.BasisUrb ;

% Candidates FACE,BY,FACE  (AFTER COMPUTING THE INTERSECTION) 
% ----------------------

[BasisINTFall,BasisRdefROT,TEXTP]= CandidatesInterfaceModesW_RVE(DATAIN,Vrb,DATAOUT,...
    nfaces,fI,BasisUdef,M,...
    TEXTP,FACES_GROUPS,SinvVal_Udef,BasisUrb,Mdom,BasisRdef,TEXTP,SinvVal_Rdef) ;

[BasisINT,TEXTP] =  ChooseIntModesWORK_RVE(BasisINTFall,SinvVal_Rdef,BasisRdefROT,TEXTP,Vrb,FACES_GROUPS,DATAIN) ; 

    DATAOUT.BasisInt = BasisINT ; 
% 
% %APPROACHES BEFORE April-21st, 2019 (before kinematical constraints)
% if DATAIN.KINEMATIC_CONSTRAINTS_MODES.ACTIVE == 0
%     [DATAOUT,V,TEXTP] =  InterfaceRVEModesStandard(DATAIN,Vrb,DATAOUT,nfaces,fI,BasisUdef,BasisRdef,M,...
%         SinvVal_Udef,SinvVal_Rdef,TEXTP,FACES_GROUPS) ;
% else
%     DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'COMBINED_ENERGETIC_METHOD',0) ;
%     if DATAIN.KINEMATIC_CONSTRAINTS_MODES.COMBINED_ENERGETIC_METHOD == 0
%         [DATAOUT,V,TEXTP] =  InterfaceRVEKinematicConstraint(DATAIN,Vrb,DATAOUT,nfaces,fI,BasisUdef,BasisRdef,M,...
%             SinvVal_Udef,SinvVal_Rdef,TEXTP,FACES_GROUPS,BasisUrb,DATA_REFMESH.M ) ;
%     else
%         % Combined with the energetic method (InterfaceRVEModesStandard)
%         DATAIN.DeformationalInterfaceModes_AlignmentMethod =1;
%         [DATAOUT,V,TEXTP] =  InterfaceRVEenergetic(DATAIN,Vrb,DATAOUT,nfaces,fI,BasisUdef,BasisRdef,M,...
%             SinvVal_Udef,SinvVal_Rdef,TEXTP,FACES_GROUPS,BasisUrb,DATA_REFMESH.M,SinvVal_Udef) ;
%         
%         
%     end
%     
% end
