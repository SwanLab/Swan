function [V_Morth,V_rbMorth] =...
    CandidatesInterfaceModes_NRIGIDB(BasisUdef_L,f1,f2,SingVal_Udef_L,...
    DATAIN,M,Vrb)

if nargin ==0
    load('tmp2.mat')
end

INCLUDE_SINGULAR_VALUES = 1; 
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
end

%BeamModesInterface = Vrb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we make the candidates for interface modes M-orthogonal to
% matrix BeamModesInterface

MAKE_M_ORTHOGONAL = 1;
if MAKE_M_ORTHOGONAL == 1
    
    Mchol = chol(M) ;
    % nmodesINTF = size(BasisUdef_L,2) ;
    % BasisREFERENCE = MorthogonalMatrix(M,Mchol,BasisREFERENCE,BeamModesInterface)  ;
    %BasisREFERENCE = BasisREFERENCE(:,1:nmodesINTF) ;  It is not advisable to
    %truncate at this level (9-Dec-2018)
    %%%%%%%%%%%%%%%
    % Now we turn BasisREFERENCE (rigid body) M-orthogonal
    Xbar = Mchol*BasisREFERENCE ;
    TOL = 0;
    [Ubar,S,Vbar] = SVDT(Xbar,TOL) ;
    V_Morth = Mchol\Ubar ;
    
    Xbar = Mchol*Vrb ;
    TOL = 0;
    [Ubar,S,Vbar] = SVDT(Xbar,TOL) ;
    V_rbMorth = Mchol\Ubar ;
    
else
     [Ubar,S,Vbar] = SVDT(BasisREFERENCE,0) ;
    V_Morth = Ubar; 
    V_rbMorth = Vrb; 
end

