function  [BasisINTFdef,BasisINTFrbORTH,BasisRdefROT,TEXTRM,BasisUdomROT,BasisINTFdef_VIRGIN  ]= ...
    CandidatesInterfaceModesSTD(DATAIN,VrbALL,DATAOUT,...
    nfaces,fI,BasisUdef,Mall,...
    TEXTP,FACES_GROUPS,SingVal_Udef_L,BasisUrb,Mdom,BasisRdef,TEXTRM,SingVal_Rdef)

if nargin == 0
    load('tmp1.mat')
end
%

DATAIN = DefaultField(DATAIN,'RotationMatrixLocal',[]) ;
hh = [] ;
LLL = [] ;
BasisINTFdef = cell(nfaces,1) ;
BasisINTFdef_VIRGIN = cell(nfaces,1) ;

BasisINTFrbORTH = cell(nfaces,1) ;

BasisRdefROT = cell(nfaces,1) ;
BasisUdomROT = cell(nfaces,1) ;

for ifgroup=1:length(FACES_GROUPS)
    f1 = fI{FACES_GROUPS{ifgroup}(1)} ;
    f2 = fI{FACES_GROUPS{ifgroup}(2)} ;
    iface1 = FACES_GROUPS{ifgroup}(1) ;
    iface2 = FACES_GROUPS{ifgroup}(2) ;
    Vrb = VrbALL{iface1} ;
    M = Mall{iface1} ;
    
    % STEP2 ) Candidate for being interface modes  (kinematically )
    % -------------------------------------
    % Vdef_Morth --> Rigid body modes (always included)
    % Vrv_Morth --> Additional modes (DEFORMATIONAL MODES)
    
    [BasisRdefROTloc,Vrb_Morth,Vdef_Morth,TEXTRM,BasisUdomROTloc] = ...
        CandidatesInterfaceModesWloc(TEXTRM,iface1,iface2,BasisRdef,f1,f2,...
        SingVal_Udef_L,BasisUdef,DATAIN,...
        Vrb,M,BasisUrb) ;
    
    BasisINTFdef_VIRGIN{iface1} = Vdef_Morth ;
    BasisINTFdef_VIRGIN{iface2} = Vdef_Morth ;
    
    % STEP 3)  Alignment method
    % --------------------------------------
    BasisRdef1 = BasisRdefROTloc{1} ;
    BasisRdef2 =  BasisRdefROTloc{2} ;
    SinvVal_Rdef  = SingVal_Rdef./SingVal_Rdef(1) ;
    
    DATAIN= DefaultField(DATAIN,'INCLUDE_SINGULAR_VALUES_COMPUTATION_OF_INTERSECTIONS',0) ; %
    if DATAIN.INCLUDE_SINGULAR_VALUES_COMPUTATION_OF_INTERSECTIONS ==1
        BasisRdef1  = bsxfun(@times,BasisRdef1',SinvVal_Rdef(1:size(BasisRdef1,2)) )' ;
        BasisRdef2   = bsxfun(@times,BasisRdef2',SinvVal_Rdef(1:size(BasisRdef1,2)))' ;
    end
    %% Orthogonal basis matrix for the intersection between BasisRdef1 and BasisRdef2  ---> R
    DATAIN= DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_AlignmentMethod',1e-4) ; %
    DATAIN= DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_ANGLES_INTERSECTION_RU',1) ; % =
    DATAIN= DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_ANGLES_INTERSECTION_V',0.001) ; % =
    
    TOL_ANGLE =DATAIN.TOL_DeformationalInterfaceModes_ANGLES_INTERSECTION_RU;
    TOL = DATAIN.TOL_DeformationalInterfaceModes_AlignmentMethod;
    TEXTRM{end+1} = '------------------------------------------------------------';
    TEXTRM{end+1} = '------------------------INTERFACE MODES-------------------------';
    TEXTRM{end+1} = ['INTERFACE =',num2str(ifgroup)];
    TEXTRM{end+1} = '------------------------------------------------------------';
    
    [R,ANGLES,sA,sB] = IntersectionSubspaces(BasisRdef1,BasisRdef2,TOL,TOL_ANGLE);
    
    TEXTRM{end+1} = ['Number of reaction modes  = ', num2str(size(BasisRdef1,2)) ];
    TEXTRM{end+1} = ['Number of reaction modes after truncation FACE -  = ', num2str(sA) ];
    TEXTRM{end+1} = ['Number of reaction modes after truncation FACE +  = ', num2str(sB) ];
    TEXTRM{end+1} = ['Dimension of reactions intersection space = ', num2str(size(R,2)) ];
    
    
    % Vdef = Vdef_Morth   ;
    %% STEP2) M-Orthogonal basis matrix for the intersection between orth. compl. of BasisUdef1 and BasisUdef2
    TEXTRM{end+1} = '------------------------------------------------------------';
    %  TEXTRM{end+1} = ['Number of def. modes  = ', num2str(size(Vdef_Morth,2)) ];
    TOL = 0 ;
    % It makes no sense to employ the concept of intersection here ....
    %[Vdef,ANGLESU,sA,sB] = IntersectionSubspaces_M(TOL,BasisUdef2,M,TOL,TOL_ANGLE) ;
    
    %%% ALIGNMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vcand= [Vrb_Morth,Vdef_Morth]  ;
    % Find h so that Vcand*h - R is minimum
    C = Vcand'*R           ;
    [U,S,V] = SVDT(C)      ;
    Vmat = Vcand*(U*V') ;
    %-------------------------------------------------
    % Intersection between  subspaces of Vdef and Vmat
    % ------------------------------------------------
    
    TOL_ANGLE =1e-3;
    TOL_SVD = 0 ;
    [Vdef_Morth,ANGLESU,sA1,sB1] = IntersectionSubspaces_M(Vdef_Morth,Vmat,M,TOL,TOL_ANGLE) ;
    
    TEXTRM{end+1} = ['Number of deformational  modes (after alignment) = ', num2str(size(Vdef_Morth,2)) ]  ;
    
    %
    
    
    
    BasisINTFdef{iface1} = Vdef_Morth ;
    BasisINTFdef{iface2} = Vdef_Morth ;
    BasisINTFrbORTH{iface1} = Vrb_Morth ;
    BasisINTFrbORTH{iface2} = Vrb_Morth ;
    BasisRdefROT{iface1} =  (BasisRdefROTloc{1}) ;
    BasisRdefROT{iface2} =  (BasisRdefROTloc{2}) ;
    
    BasisUdomROT{iface1} =  (BasisUdomROTloc{1}) ;
    BasisUdomROT{iface2} =  (BasisUdomROTloc{2}) ;
    
    
    
end


%
% figure(5623) ;
% hold on
% legend(hh,LLL)
% %
% IndexesRB = cell(nfaces,1) ;
%
% for
% end

%
% BasisRdom = cell2mat(BasisRdom) ;
% BasisUdom = cell2mat(BasisUdom) ;
