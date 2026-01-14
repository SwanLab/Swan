function [BasisREFERENCE,Vrb_Morth] =...
    CandidatesInterfaceModes_UNIFIED(BasisUdef_L,f1,f2,SingVal_Udef_L,...
    DATAIN,M,Vrb,BasisRdef,BasisUrb)

if nargin ==0
    load('tmp.mat')
end

% nfaces = 2; 
% BasisINTFrot = cell(nfaces,1) ;  % All candidates %%%%%%%%%%%
%BasisUdom = cell(nfaces,1) ;  % Original modes, but rotated 
% BasisUdom{1} = [BasisUrb(f1,:),BasisUdef_L(f1,:)] ;
% BasisUdom{2} = [BasisUrb(f2,:),BasisUdef_L(f2,:)] ;


INCLUDE_SINGULAR_VALUES = 1;  % 
BasisUdef_L1 = BasisUdef_L(f1,:) ;
BasisUdef_L2 = BasisUdef_L(f2,:) ;

SingVal_Udef_L = SingVal_Udef_L/SingVal_Udef_L(1) ;

if  INCLUDE_SINGULAR_VALUES ==1
    BasisUdef_L1 = bsxfun(@times,BasisUdef_L1',SingVal_Udef_L)' ;
    BasisUdef_L2 = bsxfun(@times,BasisUdef_L2',SingVal_Udef_L)' ;
end

DATAIN = DefaultField(DATAIN,'ROTATION_LOC',[]) ;
R = DATAIN.ROTATION_LOC ;
if isempty(R)
    BasisREFERENCE = [BasisUdef_L1, BasisUdef_L2];
else
    % Displacements of face 2 are given in domain coordinates.
    % Hence, they should be transformed back to local coordinates of the
    % interfae ---by multiplying by the transpose of the rotation matrix
    BasisUdef_L2= RotateMatrix(R',BasisUdef_L2)  ;
    BasisREFERENCE = [BasisUdef_L1, BasisUdef_L2];
  %    BasisUdom{2} =  RotateMatrix(R', BasisUdom{2})  ;
end

BeamModesInterface = Vrb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we make the candidates for interface modes M-orthogonal to
% matrix BeamModesInterface
Mchol = chol(M) ;
nmodesINTF = size(BasisUdef_L,2) ;
DATAIN = DefaultField(DATAIN,'TOLERANCE_TRUNCATE_CANDIDATE_INTERFACE_MODES',1e-6) ;  % Most reliable parameter
% for truncating the candidates. 9-Jan-2019
BasisREFERENCE = MorthogonalMatrix(M,Mchol,BasisREFERENCE,BeamModesInterface,DATAIN)  ;
%end
%%%%%%%%%%%%%%%
% %DATAIN = DefaultField(DATAIN,'SELECT_CANDIDATES_PROJECTION_REACTION_SUBSPACE',0) ; % 9-Enero-2018
% 
% if DATAIN.SELECT_CANDIDATES_PROJECTION_REACTION_SUBSPACE == 1
%     error('This option proved unreliable')
%     X =    [BasisRdef(f1,:),BasisRdef(f2,:)] ;
%     [UU,SS,VV] = SVDT(X) ;
%     % Projection onto UU
%     V_reac  = UU*(UU'*BasisREFERENCE) ;
%     [Ubar,Sbar,Vbar] =  SVDT(V_reac,0) ;
%     BasisREFERENCE = BasisREFERENCE*Vbar;
% end




% Now we turn BeamModesInterface (rigid body) M-orthogonal
Xbar = Mchol*BeamModesInterface ;
TOL = 0;
[Ubar,S,Vbar] = SVDT(Xbar,TOL) ;
Vrb_Morth = Mchol\Ubar ;



