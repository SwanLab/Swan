function [BasisINTFall,TEXTP,Indexes] =  ReactionAndInterfaceLocalEner(BasisUdef,BasisRdef,f1,f2,...
    Vrb,M,DATAIN,SinvVal_Udef,SinvVal_Rdef,iface,TEXTP,MasterDOFS_perface )
% Copy of ReactionAndInterfaceLocalModes_new  (this one was made for beam structures)
%  JAHO
if nargin == 0
    load('tmp1.mat')
end

f = [f1;f2] ;
RotationMatrixLOC = [] ;

% % Candidates to be interface fluctuation displacement modes. This function
% % also orthogonalize Vrb with respect to M
% % ----------------------------------------------
% [VrbORTH,VdefORTH] = CandidatesInterfaceModesRVE(BasisUdef,f1,f2,M,Vrb,DATAIN) ;
nmodesRB = size(Vrb,2) ;

MODES_INCLUDE =[] ;
IMODES_INCLUDE = [] ;
imode = 1;


DATAIN = DefaultField(DATAIN,'DeformationalInterfaceModes_AlignmentMethod',0) ;
DATAIN = DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_AlignmentMethod',1e-3) ;


[ BasisINT,RotationMatrixLOC,TEXTP,Vrb,Vdef,Vdef_complementary ]=  DeformModesInterface_AlignmentMethodEner(BasisRdef,f1,f2,...
    Vrb,M,DATAIN,BasisUdef,SinvVal_Udef,SinvVal_Rdef,iface,TEXTP,MasterDOFS_perface) ;

%% Matrix containing all possible modes 
% ------------------------------------------
BasisINTFall =[Vrb,Vdef,Vdef_complementary] ; 
Indexes.RB = [1:size(Vrb,2)] ;  % Rigid body DOFs
iacum = Indexes.RB(end) ;  
Indexes.DEFmaster = (iacum+1):(iacum+size(Vdef,2)) ;  % Deformational master DOFs
 iacum = iacum+size(Vdef,2) ; 

Indexes.DEFslave = (iacum+1):(iacum+size(Vdef_complementary,2)) ;   % Slave DOFs

% Rigid body modes to be included as master DOFs
DATAIN.KINEMATIC_CONSTRAINTS_MODES ...
    = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'RIGID_BODY_MODES_TO_INCLUDE',{[1:nmodesRB],[1:nmodesRB]}) ;
indRBmasterUSER = DATAIN.KINEMATIC_CONSTRAINTS_MODES.RIGID_BODY_MODES_TO_INCLUDE{iface} ; 
indRBmasterUSER = indRBmasterUSER(1:nmodesRB) ; 
% Therefore, the set of master DOFs ends up being --> 
Indexes.MASTER = [indRBmasterUSER,Indexes.DEFmaster] ; 
% And the slaves: 
rbSLAVES= setdiff(Indexes.RB,indRBmasterUSER) ; 
Indexes.SLAVES = [rbSLAVES,Indexes.DEFslave] ; 




% Rotated Reactions ---> Reactions expressed in the ref. system
% attached to the boundary interfaces, 9-Apr-2019
