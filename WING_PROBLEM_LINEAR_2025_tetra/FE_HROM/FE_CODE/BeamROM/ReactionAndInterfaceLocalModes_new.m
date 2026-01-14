function [BasisRdef,BasisINT,nBOUNDARY_INTFMODES] = ...
    ReactionAndInterfaceLocalModes_new(BasisUdef,BasisRdef,f1,f2,TOL_SINGULAR_VALUES_Hqr,...
    nBASES_BEAM,DATA_REFMESH,Vrb,M,DISP_BOUND,DATAIN,SinvVal_Udef,SinvVal_Rdef)
if nargin == 0
    load('tmp4.mat')
end
 
f = [f1;f2] ;
% We have distinguish between "beam" modes and remaining modes (local effects)
%
BasisUdef_B = BasisUdef(:,1:nBASES_BEAM.DISPLACEMENTS) ; % Displacements, beam modes

BasisUdef_L = BasisUdef(:,nBASES_BEAM.DISPLACEMENTS+1:end) ; % Displacements, local effects
SingVal_Udef_L = SinvVal_Udef(nBASES_BEAM.DISPLACEMENTS+1:end)  ; % Associated singular values
BasisRdef_L = BasisRdef(:,nBASES_BEAM.REACTIONS+1:end) ;     % Reactions, local effects
if ~isempty(SinvVal_Rdef)
SingVal_Rdef_L = SinvVal_Rdef(nBASES_BEAM.REACTIONS+1:end)  ; % Associated singular values
else
    SingVal_Rdef_L = 1 ; 
end

% What to do with these modes  ?
%  Reference matrix

% Candidates to be interface displacement modes
% ----------------------------------------------
if ~isempty(DATA_REFMESH.RotationMatrixFace)
iface = 2; 
DATAIN.ROTATION_LOC = DATA_REFMESH.RotationMatrixFace{iface} ; 
R = DATAIN.ROTATION_LOC ; 
else
    R = [] ; 
end


[nBOUNDARY_INTFMODES,BasisREFERENCE,Vrb_Morth] = CandidatesInterfaceModes(BasisUdef_L,f1,f2,DISP_BOUND,SingVal_Udef_L,...
    DATAIN,BasisUdef_B,M,Vrb) ;

% Determining number of reaction modes
DATAIN = DefaultField(DATAIN,'MethodForSelectingReactionModes',1) ;

if  DATAIN.MethodForSelectingReactionModes == 0
    % Method used before 29-09-2018
    BasisRdef_L = ReactionModesLocalEffects(BasisUdef_L,BasisRdef_L,f,TOL_SINGULAR_VALUES_Hqr) ;
else
    BasisRdef_L = ReactionModesLocalEffectsNEW(BasisUdef_L,BasisRdef_L,f,TOL_SINGULAR_VALUES_Hqr,SingVal_Rdef_L) ;
    
end

BasisRdef = [BasisRdef(:,1:nBASES_BEAM.REACTIONS), BasisRdef_L] ;

DATAOUT.BasisRdef = BasisRdef ;
% --------
%% Final selection of interface modes 
% -------------------------------------
OLD_METHOD  = 0; % Method that only compares local effects modes (version previous to December 9th-2018)

if OLD_METHOD == 1
BasisINT = SelectCandidatesInterfaceModes_BEAMS(Vrb,M,BasisRdef_L,BasisREFERENCE,...
    TOL_SINGULAR_VALUES_Hqr,f1,f2,R,DATAIN) ; 
else
    BasisRdef_L = BasisRdef ; 
   BasisINT = SelectCandidatesInterfaceModes_BEAMSnew(Vrb,M,BasisRdef_L,BasisREFERENCE,...
    TOL_SINGULAR_VALUES_Hqr,f1,f2,R,DATAIN,Vrb_Morth) ;  
end