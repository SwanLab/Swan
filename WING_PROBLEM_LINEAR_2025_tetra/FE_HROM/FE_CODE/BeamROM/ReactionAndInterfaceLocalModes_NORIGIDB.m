function [BasisRdef,BasisINT] = ...
    ReactionAndInterfaceLocalModes_NORIGIDB(BasisUdef,BasisRdef,f1,f2,TOL_SINGULAR_VALUES_Hqr,...
    nBASES_BEAM,DATA_REFMESH,Vrb,M,DATAIN,SingVal_Udef,SingVal_Rdef)
if nargin == 0
    load('tmp.mat')
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

% Candidate for being interface modes
% Vdef_Morth --> Rigid body modes (always included)
% Vrv_Morth --> Additional modes
[V_Morth,Vrb_Morth] = CandidatesInterfaceModes_NRIGIDB(BasisUdef,f1,f2,SingVal_Udef,...
    DATAIN,M,Vrb) ;

 % Selection of reaction modes 
 % ---------------------------
BasisRdef = ReactionModesLocalEffectsNEW(BasisUdef,BasisRdef,f,TOL_SINGULAR_VALUES_Hqr,SingVal_Rdef,size(Vrb,2)) ;
  
DATAOUT.BasisRdef = BasisRdef ;
% --------
%% Final selection of interface modes 
BasisINT = SelectCandidatesInterfaceModes_NRBODY(Vrb,M,BasisRdef,V_Morth,...
    TOL_SINGULAR_VALUES_Hqr,f1,f2,R,DATAIN,Vrb_Morth) ;
 