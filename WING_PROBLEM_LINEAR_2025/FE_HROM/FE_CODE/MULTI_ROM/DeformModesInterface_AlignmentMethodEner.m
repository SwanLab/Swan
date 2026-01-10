function [BasisINT,RotationMatrixLOC,TEXTP,Vrb,Vdef,Vdef_complementary ]...
    =  DeformModesInterface_AlignmentMethodEner(BasisRdef,f1,f2,...
    Vrb,M,DATAIN,BasisUdef,SinvVal_Udef,SinvVal_Rdef,iface,TEXTP,...
    MasterDOFS_perface)
if nargin == 0
    load('tmp.mat')
end

if size(Vrb,2) == 3
    ndim = 2 ;
else
    ndim = 3 ;
end





% nmodesRB = size(Vrb,2) ;
% Vcandidate = [VrbORTH,VdefORTH] ;

% STEP 1: Rotation of matrices (def. and reactions. )
% --------------------------------------------------------------------------------------------
DATAIN = DefaultField(DATAIN,'RotationMatrixLocal',[]) ;

if isempty(DATAIN.RotationMatrixLocal)
    %     ReactionsALL = [BasisRdef(f1,:),BasisRdef(f2,:)] ; % Reactions of each pair of bound. interfaces
    %     DispALL = [BasisUdef(f1,:),BasisUdef(f2,:)] ;
    BasisRdef1 = BasisRdef(f1,:) ;
    BasisUdef1 = BasisUdef(f1,:) ;
    
    BasisRdef2= BasisRdef(f2,:) ;
    BasisUdef2 = BasisUdef(f2,:) ;
    RotationMatrixLOC = {eye(ndim),eye(ndim)} ;
    
else
    R1 = DATAIN.RotationMatrixLocal{1} ;
    R2 = DATAIN.RotationMatrixLocal{2} ;
    RotationMatrixLOC = cell(1,2) ;
    if ~isempty(R1)
        BasisRdef1= RotateMatrix(R1,BasisRdef(f1,:))  ;
        BasisUdef1= RotateMatrix(R1,BasisUdef(f1,:))  ;
        RotationMatrixLOC{1} = R1;
    else
        BasisRdef1 = BasisRdef(f1,:) ;
        BasisUdef1 = BasisUdef(f1,:) ;
        RotationMatrixLOC{1} = eye(ndim);
    end
    
    if ~isempty(R2)
        BasisRdef2= RotateMatrix(R2,BasisRdef(f2,:))  ;
        BasisUdef2= RotateMatrix(R2,BasisUdef(f2,:))  ;
        RotationMatrixLOC{2} = R2;
    else
        BasisRdef2 = BasisRdef(f2,:) ;
        BasisUdef2 = BasisUdef(f2,:) ;
        RotationMatrixLOC{2} = eye(ndim);
    end
    % ReactionsALL = [BasisRdef1,BasisRdef2] ;
    %   ReactionsALL = [BasisRdef1,BasisRdef2] ;
end

%RotatedReactionsT = {BasisRdef1',BasisRdef2'} ;

SinvVal_Rdef  = SinvVal_Rdef./SinvVal_Rdef(1) ;
SinvVal_Udef = SinvVal_Udef./SinvVal_Udef(1) ;

DATAIN= DefaultField(DATAIN,'INCLUDE_SINGULAR_VALUES_COMPUTATION_OF_INTERSECTIONS',0) ; %
if DATAIN.INCLUDE_SINGULAR_VALUES_COMPUTATION_OF_INTERSECTIONS ==1
    BasisRdef1  = bsxfun(@times,BasisRdef1',SinvVal_Rdef(1:size(BasisRdef1,2)) )' ;
    BasisRdef2   = bsxfun(@times,BasisRdef2',SinvVal_Rdef(1:size(BasisRdef1,2)))' ;
%       BasisUdef1 = bsxfun(@times,BasisUdef1',SinvVal_Udef(1:size(BasisUdef1,2)))' ;
%      BasisUdef2= bsxfun(@times,BasisUdef2',SinvVal_Udef(1:size(BasisUdef1,2)))' ;
end

%% Orthogonal basis matrix for the intersection between BasisRdef1 and BasisRdef2  ---> R
DATAIN= DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_AlignmentMethod',1e-4) ; %
DATAIN= DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_ANGLES_INTERSECTION_RU',1) ; % =
DATAIN= DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_ANGLES_INTERSECTION_V',0.001) ; % =
DATAIN= DefaultField(DATAIN,'USE_INTERSECTIONS_FOR_CANDIDATES',1) ; % =


 
TOL_ANGLE =DATAIN.TOL_DeformationalInterfaceModes_ANGLES_INTERSECTION_RU;
TOL = DATAIN.TOL_DeformationalInterfaceModes_AlignmentMethod;
TEXTP{end+1} = '------------------------------------------------------------';
TEXTP{end+1} = '------------------------INTERFACE MODES-------------------------';
TEXTP{end+1} = ['INTERFACE =',num2str(iface)];
TEXTP{end+1} = '------------------------------------------------------------';
TEXTP{end+1} = ['Number of reaction modes  = ', num2str(size(BasisRdef1,2)) ];
DATALOC.RELATIVE_SVD = 1 ;


if DATAIN.USE_INTERSECTIONS_FOR_CANDIDATES == 1
    [QA,SA,~] = SVDT(BasisRdef1,TOL,DATALOC) ;
    [QB,SB,~] = SVDT(BasisRdef2,TOL,DATALOC) ;
    TEXTP{end+1} = ['Number of reaction modes after truncation FACE -  = ', num2str(size(QA,2)) ];
    TEXTP{end+1} = ['Number of reaction modes after truncation FACE +  = ', num2str(size(QB,2)) ];
    [Y,CT,Z] = SVDT(QA'*QB) ;
    ANGLES = real(acos(CT))*180/pi ;
    nREACT = length(find(ANGLES <= TOL_ANGLE));
    R = QA*Y(:,1:nREACT) ;
    TEXTP{end+1} = ['Dimension of reactions intersection space = ', num2str(size(R,2)) ];
else
    [R,SA,~] = SVDT([BasisRdef1,BasisRdef2],TOL,DATALOC) ;
    TEXTP{end+1} = ['Number of reactions = ', num2str(size(R,2)) ];
end



%RotatedReactionsT = {R',R'} ;
% Replace it  by ---> R = IntersectionSubspaces(A,B,TOL_SVD,TOL_ANGLE)
% !!!

%% STEP2) M-Orthogonal basis matrix for the intersection between orth. compl. of BasisUdef1 and BasisUdef2
TEXTP{end+1} = '------------------------------------------------------------';
TEXTP{end+1} = ['Number of displacement modes  = ', num2str(size(BasisUdef1,2)) ];

PG = (Vrb'*M*Vrb)  ;
BasisUdef1 = BasisUdef1  - Vrb*(PG\(Vrb'*M*BasisUdef1)) ;
BasisUdef2 = BasisUdef2  - Vrb*(PG\(Vrb'*M*BasisUdef2)) ;

BasisUdefALL_orthM = [BasisUdef1,BasisUdef2]  ;  % Orthogonal to M 
Mchol = chol(M) ;

if DATAIN.USE_INTERSECTIONS_FOR_CANDIDATES == 1
 
BasisUdef1M = Mchol*BasisUdef1 ;
BasisUdef2M = Mchol*BasisUdef2 ;
[QA,SA,~] = SVDT(BasisUdef1M,TOL,DATALOC) ;
[QB,SB,~] = SVDT(BasisUdef2M,TOL,DATALOC) ;
TEXTP{end+1} = ['Number of disp. modes after truncation FACE -  = ', num2str(size(QA,2)) ];
TEXTP{end+1} = ['Number of disp modes after truncation FACE +  = ', num2str(size(QB,2)) ];


[Y,CT,Z] = SVDT(QA'*QB) ;
ANGLES = real(acos(CT))*180/pi ;
%TOL_ANGLE = 1;
nDISP = length(find(ANGLES <= TOL_ANGLE));
Vdef = Mchol\(QA*Y(:,1:nDISP)) ;

TEXTP{end+1} = ['Dimension of displacem. intersection space = ', num2str(size(Vdef,2)) ];


else
    DATALOC.TOL = TOL ; 
  [Vdef,SA,~] = WSVDT(BasisUdefALL_orthM,M,DATALOC) ;
  
end

% WORK DONE BY THE CANDIDATES Vdef 

% FACE 1 
Vdef_norm_inv = 1./sqrt(sum(Vdef.^2,1)) ;
Vdef_unit = bsxfun(@times,Vdef',Vdef_norm_inv')' ; 

BasisRdef1_norm_inv = 1./sqrt(sum(BasisRdef1.^2,1)) ;
BasisRdef1_unit = bsxfun(@times,BasisRdef1',BasisRdef1_norm_inv')' ; 

BasisRdef2_norm_inv = 1./sqrt(sum(BasisRdef2.^2,1)) ;
BasisRdef2_unit = bsxfun(@times,BasisRdef2',BasisRdef2_norm_inv')' ; 

W1 = Vdef_unit'*BasisRdef1_unit ; 
W2 = Vdef_unit'*BasisRdef2_unit ; 


 

% REplace it by R = IntersectionSubspaces_M(A,B,M,TOL_SVD,TOL_ANGLE)
% !!!!!!!!!

%%%%%

% % Orthogonal complement
% PG = (Vrb'*M*Vrb)  ;
% Rast = R  - Vrb*(PG\(Vrb'*M*R)) ;

VrbORTH = WSVDT(Vrb,M) ;
Vcand= [VrbORTH,Vdef]  ;
% Find h so that Vcand*h - R is minimum
C = Vcand'*R           ;
[U,S,V] = SVDT(C)      ;

Vmat = Vcand*(U*V') ;

%-------------------------------------------------
% Intersection between  subspaces of Vdef and Vmat
% ------------------------------------------------
[QA,~,~ ]= SVDT(Mchol*Vdef)    ;
[QB,~,~ ]= SVDT(Mchol*Vmat)    ;
C = QA'*QB                     ;
[Y,CTHETA,Z] = SVDT(C)         ;
% % CTHETA --> Cosine principal angles
TOL_ANGLE = DATAIN.TOL_DeformationalInterfaceModes_ANGLES_INTERSECTION_V ;
ANGLES = real(acos(CTHETA))*180/pi      ;
nDEF = length(find(ANGLES<= TOL_ANGLE)) ;

%nDEF =min(nDEF,MasterDOFS_perface) ;  % MasterDOFS_perface ---> By default, a very high value

TEXTP{end+1} = ['Number of deformational displacements (computed) = ', num2str(nDEF) ];


DATAIN = DefaultField(DATAIN,'NUMBER_DEFORMATIONAL_MODES_FACES_13_24',[])  ;
if ~isempty(DATAIN.NUMBER_DEFORMATIONAL_MODES_FACES_13_24)
    nDEFpres  =DATAIN.NUMBER_DEFORMATIONAL_MODES_FACES_13_24(iface) ;
    if ~isempty(nDEFpres)
        nDEF = min(nDEFpres,nDEF) ;
    end
end
TEXTP{end+1} = ['Number of deformational displacements (final, user) = ', num2str(nDEF) ];


Vdef = Mchol\QA*Y(:,1:nDEF)             ;
BasisINT = [Vrb,Vdef]                   ;

%%%%%%%%%%%%%%%%%%%%%%%%%
if  ~isempty(Vdef)
    % Orthogonal complement
    PG = (Vdef'*M*Vdef)  ;
    BasisUdefALL_orthM = BasisUdefALL_orthM  - Vdef*(PG\(Vdef'*M*BasisUdefALL_orthM)) ;
    
end
DATALOC.TOL = 1e-6 ;
[Vdef_complementary,S,~] = WSVDT(BasisUdefALL_orthM,M,DATALOC) ;

% Therefore, we have that [BasisINT,Vdef_complementary]  span all possible
% configurations of these faces. So, in the terminology of the paper, this
% means that, for face (-), the master DOFs should be among
% 1) The rigid body modes  (1,2, ...nRB)
% 2) The energetically dominant deformational modes (nRB+1,nRB+2,...nRB+nDEF)
% On the other hand, the slave DOFs would be formed by those not selected
% among these, plus the columns nRB+nDEF+1,nRB+nDEF+2 .... 




% BasisINTprueba = [VrbORTH,Vdef] ;
% BasisINT = BasisINTprueba ;

%    [UVer1,Sver1,Vver1 ] = SVDT(BasisINTprueba'*BasisRdef1) ;
%    [UVer2,Sver2,Vver2 ] = SVDT(BasisINTprueba'*BasisRdef2) ;
