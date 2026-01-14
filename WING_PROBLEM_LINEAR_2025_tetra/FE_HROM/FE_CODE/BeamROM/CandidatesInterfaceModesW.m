function [BasisREFERENCE,Vrb_Morth,TEXTRM,BasisRdefROT] =...
    CandidatesInterfaceModesW(BasisUdef_L,f1,f2,SingVal_Udef_L,...
    DATAIN,M,Vrb,TEXTRM,BasisRdef)

if nargin ==0
    load('tmp1.mat')
end
 
BasisRdefROT = {BasisRdef(f1,:); BasisRdef(f2,:)} ; 

BasisUdef_L1 = BasisUdef_L(f1,:) ;
BasisUdef_L2 = BasisUdef_L(f2,:) ;

SingVal_Udef_L = SingVal_Udef_L/SingVal_Udef_L(1) ;

BasisUdef_L1 = bsxfun(@times,BasisUdef_L1',SingVal_Udef_L)' ;
BasisUdef_L2 = bsxfun(@times,BasisUdef_L2',SingVal_Udef_L)' ;
 
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
    BasisRdefROT{2} = RotateMatrix(R',  BasisRdefROT{2} )  ;
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
 

TEXTRM{end+1} = ['********************************************'] ; 
TEXTRM{end+1} = ['Total number of candidates for INTEFACE MODES =  ',num2str(size(Vrb,2)),' + ',num2str(size(BasisREFERENCE,2)),' = ',...
    num2str(size(Vrb,2) + size(BasisREFERENCE,2))] ; 

% Now we turn BeamModesInterface (rigid body) M-orthogonal
Xbar = Mchol*BeamModesInterface ;
TOL = 0;
[Ubar,S,Vbar] = SVDT(Xbar,TOL) ;
Vrb_Morth = Mchol\Ubar ;



