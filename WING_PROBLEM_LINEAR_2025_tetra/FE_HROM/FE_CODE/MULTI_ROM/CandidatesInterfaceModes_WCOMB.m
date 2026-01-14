function  [BasisINTFdefCAND,BasisRdom,BasisUdom,TEXTRM,BasisINTFrbORTH]= CandidatesInterfaceModes_WCOMB(DATAIN,VrbALL,DATAOUT,...
    nfaces,fI,BasisUdef,Mall,...
    TEXTRM,FACES_GROUPS,SingVal_Udef_L,BasisUrb,Mdom,BasisRdef,SingVal_Rdef)

if nargin == 0
    load('tmp1.mat')
end

BasisINTFdefCAND = cell(nfaces,1) ;
BasisINTFrbORTH = cell(nfaces,1) ;

BasisRdom = cell(nfaces,1) ;
BasisUdom = cell(nfaces,1) ;

DATALOC = DATAIN.KINEMATIC_CONSTRAINTS_MODES  ; 

DATAIN = DefaultField(DATAIN,'RotationMatrixLocal',[]) ;
% Tolerance determining intersection pair of faces
%DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'COMPUTE_INTERSECTION',1) ;
%COMPUTE_INTERSECTION =DATAIN.KINEMATIC_CONSTRAINTS_MODES.COMPUTE_INTERSECTION ;
%DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'TOL_INTERSECT_DISP',0.1) ;
DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'TOL_SVD_DISP_MODES',1e-4) ;

for ifgroup=1:length(FACES_GROUPS)
    f1 = fI{FACES_GROUPS{ifgroup}(1)} ;
    f2 = fI{FACES_GROUPS{ifgroup}(2)} ;
    iface1 = FACES_GROUPS{ifgroup}(1) ;
    iface2 = FACES_GROUPS{ifgroup}(2) ;
    Vrb = VrbALL{iface1} ;
    M = Mall{iface1} ;
    
    TEXTRM{end+1} = ['--------------------------------------------'] ;
    TEXTRM{end} = ['CANDIDATES FOR BEING INTERFACE MODES']  ;
    TEXTRM{end+1} = ['faces = ',num2str(iface1),' and ',num2str(iface2)]  ;
    
    % Reactive forces (resultants and selfequilibrated)
    BasisRrb1 = Mdom(f1,f1)*BasisUrb(f1,:) ;
    BasisRrb2 = Mdom(f2,f2)*BasisUrb(f2,:) ;
    BasisRdef1 =  BasisRdef(f1,:) ;
    BasisRdef2 =  BasisRdef(f2,:) ;
    
    BasisRdom1 = [BasisRrb1,BasisRdef1] ;
    BasisRdom2 = [BasisRrb2,BasisRdef2] ;
    
    BasisUdom1 = [BasisUrb(f1,:),BasisUdef(f1,:)] ;
    BasisUdom2 = [BasisUrb(f2,:),BasisUdef(f2,:)] ;
    
    if isempty(DATAIN.RotationMatrixLocal)
        
        BasisUdef1 = BasisUdef(f1,:) ;
        BasisUdef2 = BasisUdef(f2,:) ;
        
        
    else
        error('Option not implemented')
        
        R1 = DATAIN.RotationMatrixLocal{1} ;
        R2 = DATAIN.RotationMatrixLocal{2} ;
        if ~isempty(R1)
            BasisUdef1= RotateMatrix(R1',BasisUdef(f1,:))  ;
            BasisUdef2= RotateMatrix(R2',BasisUdef(f2,:))  ;
            BasisRdom1= RotateMatrix(R1',BasisRdom1)  ;
            BasisRdom2= RotateMatrix(R2',BasisRdom2)  ;
            
            BasisUdom1= RotateMatrix(R1',BasisUdom1)  ;
            BasisUdom2= RotateMatrix(R2',BasisUdom2)  ;
        else
            BasisUdef1 = BasisUdef(f1,:) ;
            BasisUdef2 = BasisUdef(f2,:) ;
        end
    end
    
    
    [BasisRdefROTloc,Vrb_Morth,Vdef_Morth,TEXTRM] = ...
        CandidatesInterfaceModesWloc(TEXTRM,iface1,iface2,BasisRdef,f1,f2,SingVal_Udef_L,BasisUdef,DATAIN,...
        Vrb,M,BasisUrb) ;
    
    % STEP 3)  INTERSECTION OF SUBSPACES
    % --------------------------------------
    DATALOC = DefaultField(DATALOC,'TOLERANCE_SVD_REACTION_MODES_INTERSECTION',1e-6) ; 
    DATALOC = DefaultField(DATALOC,'TOLERANCE_ANGLE_INTERSECTION_REACTIONS',0.1) ; 
    TOL_SVD = DATALOC.TOLERANCE_SVD_REACTION_MODES_INTERSECTION;
    TOL_ANGLE = DATALOC.TOLERANCE_ANGLE_INTERSECTION_REACTIONS ;
    [RR,ANGLES ]= IntersectionSubspaces(BasisRdefROTloc{1},BasisRdefROTloc{2},TOL_SVD,TOL_ANGLE) ;
    nmodesR = size(RR,2) ;
    
%     figure(5623) ;
%     hold on
%     hh(end+1) = plot(ANGLES) ;
%     xlabel('Modes') ;
%     ylabel('Angles (degrees)') ;
%     title(['INTERSECTION BETWEEN SUBSPACES BasisRdef_{f1} and BasisRdef_{f2}'])  ;
%     LLL{end+1} = ['FACES ',num2str(iface1),'  ',num2str(iface2)] ;
%     TEXTRM{end+1} = ['Dimension of the intersection of reaction subspaces (for TOL ANGLE = ',num2str(TOL_ANGLE),') --> ',num2str(nmodesR)] ;
%     TEXTRM{end+1} = '****************************************' ;
    % STEP 4)  Selection of candidates
    % ----------------------------------
    [BasisINT,TEXTRM] = ...
        ChooseIntModesWORK1_RVE(Vrb_Morth,Vdef_Morth,SingVal_Rdef,BasisRdefROTloc,TEXTRM,nmodesR,Vrb) ;
    
    nRB = size(Vrb,2) ; 
    Vdef = BasisINT(:,nRB+1:end) ; 
    
     
    
    %
%     %
%     %% STEP2) M-Orthogonal basis matrix for the intersection between orth. compl. of BasisUdef1 and BasisUdef2
%     TEXTP{end+1} = '------------------------------------------------------------';
%     TEXTP{end+1} = ['Number of displacement modes  = ', num2str(size(BasisUdef1,2)) ];
%     
%     PG = (Vrb'*M*Vrb)  ;
%     BasisUdef1 = BasisUdef1  - Vrb*(PG\(Vrb'*M*BasisUdef1)) ;
%     BasisUdef2 = BasisUdef2  - Vrb*(PG\(Vrb'*M*BasisUdef2)) ;
    
%     
%     % if COMPUTE_INTERSECTION == 0
%     % No intersection. Just a M-orthtogonal basis matrix
%     %        error('This option should be properly examined. It produces ill-conditioned equations')
%     DATALOC.RELATIVE_SVD = 1 ;
%     DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'INCLUDE_SINGULAR_VALUES_SVD_NOINTER',1) ;
%     if DATAIN.KINEMATIC_CONSTRAINTS_MODES.INCLUDE_SINGULAR_VALUES_SVD_NOINTER == 1
%         BasisUdef1_SV = bsxfun(@times,BasisUdef1',SinvVal_Udef)' ;
%         BasisUdef2_SV = bsxfun(@times,BasisUdef2',SinvVal_Udef)' ;
%     else
%         BasisUdef1_SV = BasisUdef1; BasisUdef2_SV = BasisUdef2 ;
%     end
%     DATALOC.TOL =DATAIN.KINEMATIC_CONSTRAINTS_MODES.TOL_SVD_DISP_MODES ;;
%     [Vdef,SV,~] = WSVDT([BasisUdef1_SV BasisUdef2_SV],M,DATALOC) ;
%     TEXTP{end+1} = ['Number of def. modes   = ', num2str(size(BasisUdef1,2)) ];
    
    
    
    
    %   else
    %         Mchol = chol(M) ;
    %         BasisUdef1M = Mchol*BasisUdef1 ;
    %         BasisUdef2M = Mchol*BasisUdef2 ;
    %         DATALOC.RELATIVE_SVD = 1 ;
    %         TOL = DATAIN.KINEMATIC_CONSTRAINTS_MODES.TOL_SVD_DISP_MODES ;
    %         [QA,SA,~] = SVDT(BasisUdef1M,TOL,DATALOC) ;
    %         [QB,SB,~] = SVDT(BasisUdef2M,TOL,DATALOC) ;
    %         TEXTP{end+1} = ['Number of disp. modes after truncation FACE -  = ', num2str(size(QA,2)) ];
    %         TEXTP{end+1} = ['Number of disp modes after truncation FACE +  = ', num2str(size(QB,2)) ];
    %
    %         [Y,CT,Z] = SVDT(QA'*QB) ;
    %         ANGLES = real(acos(CT))*180/pi ;
    %         TOL_ANGLE = DATAIN.KINEMATIC_CONSTRAINTS_MODES.TOL_INTERSECT_DISP; %,0.1) ;
    %         nDISP = length(find(ANGLES <= TOL_ANGLE));
    %         Vdef = Mchol\(QA*Y(:,1:nDISP)) ;
    
   % TEXTP{end+1} = ['Dimension of displacem. intersection space = ', num2str(size(Vdef,2)) ];
    %  end
    
    
    BasisINTFdefCAND{iface1} = Vdef ;
    BasisINTFdefCAND{iface2} = Vdef ;
    BasisINTFrbORTH{iface1} = Vrb_Morth ; 
    BasisINTFrbORTH{iface2} = Vrb_Morth ; 
    BasisRdom{iface1} =  (BasisRdom1) ;
    BasisRdom{iface2} =  (BasisRdom2) ;
    BasisUdom{iface1} =  (BasisUdom1) ;
    BasisUdom{iface2} =  (BasisUdom2) ;
end

BasisRdom = cell2mat(BasisRdom) ;
BasisUdom = cell2mat(BasisUdom) ;
