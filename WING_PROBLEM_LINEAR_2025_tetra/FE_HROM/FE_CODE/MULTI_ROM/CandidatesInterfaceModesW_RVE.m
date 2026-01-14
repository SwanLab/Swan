function  [BasisINTFall,BasisRdefROT,TEXTRM  ]= CandidatesInterfaceModesW_RVE(DATAIN,VrbALL,DATAOUT,...
    nfaces,fI,BasisUdef,Mall,...
    TEXTP,FACES_GROUPS,SingVal_Udef_L,BasisUrb,Mdom,BasisRdef,TEXTRM,SingVal_Rdef)

if nargin == 0
    load('tmp1.mat')
end
%
DATALOC = DATAIN.INTERFACE_MODES_REACTIONS_WORK ;

DATALOC = DefaultField(DATALOC,'TOLERANCE_SVD_DEF_MODES',1e-6) ; % TOLERANCE FOR SVD DEFORMATIONAL MODES
DATALOC = DefaultField(DATALOC,'TOLERANCE_ANGLE_INTERSECTION_REACTIONS',0.01) ; % TOLERANCE FOR SVD DEFORMATIONAL MODES
DATALOC = DefaultField(DATALOC,'TOLERANCE_SVD_REACTION_MODES_INTERSECTION',1e-6) ; % TOLERANCE FOR SVD reaction MODES

DATAIN = DefaultField(DATAIN,'RotationMatrixLocal',[]) ;
hh = [] ;
LLL = [] ;
BasisINTFall = cell(nfaces,1) ;
BasisRdefROT = cell(nfaces,1) ;

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
    
    [BasisRdefROTloc,Vrb_Morth,Vdef_Morth,TEXTRM] = ...
        CandidatesInterfaceModesWloc(TEXTRM,iface1,iface2,BasisRdef,f1,f2,SingVal_Udef_L,BasisUdef,DATAIN,...
        Vrb,M) ;
    
    % STEP 3)  INTERSECTION OF SUBSPACES
    % --------------------------------------
    TOL_SVD = DATALOC.TOLERANCE_SVD_REACTION_MODES_INTERSECTION;
    TOL_ANGLE = DATALOC.TOLERANCE_ANGLE_INTERSECTION_REACTIONS ;
    [RR,ANGLES ]= IntersectionSubspaces(BasisRdefROTloc{1},BasisRdefROTloc{2},TOL_SVD,TOL_ANGLE) ;
    nmodesR = size(RR,2) ;
    
    figure(5623) ;
    hold on
    hh(end+1) = plot(ANGLES) ;
    xlabel('Modes') ;
    ylabel('Angles (degrees)') ;
    title(['INTERSECTION BETWEEN SUBSPACES BasisRdef_{f1} and BasisRdef_{f2}'])  ;
    LLL{end+1} = ['FACES ',num2str(iface1),'  ',num2str(iface2)] ;
    TEXTRM{end+1} = ['Dimension of the intersection of reaction subspaces (for TOL ANGLE = ',num2str(TOL_ANGLE),') --> ',num2str(nmodesR)] ;
    TEXTRM{end+1} = '****************************************' ;
    % STEP 4)  Selection of candidates
    % ----------------------------------
    [BasisINT,TEXTRM] = ...
        ChooseIntModesWORK1_RVE(Vrb_Morth,Vdef_Morth,SingVal_Rdef,BasisRdefROTloc,TEXTRM,nmodesR,Vrb) ;
    
    
    
    
    BasisINTFall{iface1} = BasisINT ;
    BasisINTFall{iface2} = BasisINT ;
    BasisRdefROT{iface1} =  (BasisRdefROTloc{1}) ;
    BasisRdefROT{iface2} =  (BasisRdefROTloc{2}) ;
    
end

figure(5623) ;
hold on
legend(hh,LLL)
% 
% IndexesRB = cell(nfaces,1) ; 
%  
% for  
% end

%
% BasisRdom = cell2mat(BasisRdom) ;
% BasisUdom = cell2mat(BasisUdom) ;
