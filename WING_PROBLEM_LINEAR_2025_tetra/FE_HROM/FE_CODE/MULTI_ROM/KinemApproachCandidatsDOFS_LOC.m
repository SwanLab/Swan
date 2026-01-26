function [TEXTP,DATAOUT,V]  =KinemApproachCandidatsDOFS_LOC(BasisINTFdefCAND,VrbALL,fI,...
    BasisUrb,BasisUdef,Mdom,TEXTP,DATAOUT,BasisINTFrbORTH,BasisUdomROT,DATAIN,Minterfaces)

if nargin == 0
    load('tmp2.mat')
   % DATAIN.KinematicCONSTR_OnlyTranslatRB_DOFS = 1;
    DATAIN.MATRIX_TO_USE_IN_MINIMIZATION_DMATRIX = 'ROTATED_INTERFACES' ; 
end

[~,nDEF] = cellfun(@size,BasisINTFdefCAND) ;  %Number of deformational modes for each interface
[~,nRB] = cellfun(@size,VrbALL) ;   % Number of RB modes for each interface

nCAND = nDEF(:) + nRB(:) ;  % Total number of candidates for each interface

% STEP 2. All candidates
% ----------------------
BasisINTFcand =  cell(size(BasisINTFdefCAND)) ;
BasisINTFcandORTH = cell(size(BasisINTFdefCAND)) ;
IndicesRB = [] ;
%IndicesDEF = [] ;
iacum = 0 ;
% MAtrix with all the candidates, including RB modes
for iface = 1:length(VrbALL)
    IndicesRBloc =  (1:nRB(iface))+iacum ;
    
    IndicesRB = [IndicesRB,IndicesRBloc] ;
    BasisINTFcand{iface} = sparse([VrbALL{iface},BasisINTFdefCAND{iface}]) ;
    BasisINTFcandORTH{iface} = [BasisINTFrbORTH{iface},BasisINTFdefCAND{iface}];
    %     IndicesDEFloc =  (nRB(iface):size(BasisINTFcand{iface},2)  )+iacum ;
    %      IndicesDEF{end+1} =
    iacum = iacum +  nRB(iface) + nDEF(iface) ;
end
[nelems,nDOFsFACEall] = cellfun(@size,BasisINTFcand) ;

IndicesDEF = 1:sum(nCAND) ;
IndicesDEF(IndicesRB) = [] ;

% STEP 3. Form a diagonal basis matrix  (this is matrix U in the paper)
% ------------------------------------
BasisINTFall = blkdiag(BasisINTFcand{:}) ;

% STEP 4. Basis of domain displacements
f = cell2mat(fI') ;
% BasisUrb_f =BasisUrb(f,:) ;
% BasisUdef_f = BasisUdef(f,:) ;


% BasisUdom_f = [BasisUrb_f,BasisUdef_f] ;
%
% BasisUDEFf12 =[  BasisUdef([fI{1}; fI{2}],:)];


% STEP 5. Solving minimization problem



% This is equivalent to solve a least-squares problem in the norm of
% the mass matrix. Pertinent tests were made and the results appeared
% to confirm that the implementation is in principle correct

  % Shouldn't this matrix be expressed in the reference system of the
% interfaces ?
DATAIN= DefaultField(DATAIN,'MATRIX_TO_USE_IN_MINIMIZATION_DMATRIX','ROTATED_INTERFACES') ; % = 'DOMAIN' ; % 'ROTATED_INTERFACES'
switch DATAIN.MATRIX_TO_USE_IN_MINIMIZATION_DMATRIX
    case 'DOMAIN'
        Mchol =   chol(Mdom(f,f)) ;
        % Implementation before November 15-th 2019
    case 'ROTATED_INTERFACES'
        M = blkdiag(Minterfaces{:}) ;
        Mchol = chol(M) ; 
        
end



MBasisINTFall = Mchol*BasisINTFall ;
% MBasisUdom = Mchol*BasisUdom_f ;
BasisUdomROT = cell2mat(BasisUdomROT) ;
MBasisUdom = Mchol*BasisUdomROT ;


Dcomp = MBasisINTFall\MBasisUdom ;

% IndicesRB_input = IndicesRB ;
% DATAIN = DefaultField(DATAIN,'KinematicCONSTR_OnlyTranslatRB_DOFS',0) ;
% if DATAIN.KinematicCONSTR_OnlyTranslatRB_DOFS == 1;
%     nfaces = length(nRB) ;
%     IndicesRB_input = reshape(IndicesRB,[],nfaces) ;
%     if nRB(1) == 3
%         % 2D problems
%         TRANSIND = [1:2] ;
%     else
%         TRANSIND = [1:3] ;
%     end
%     IndicesRB_input = IndicesRB_input(TRANSIND,:) ;
%     IndicesRB_input = IndicesRB_input(:)' ;
% end

[DOFm,DOFs] = CoarseDOFS_DEIMbased(Dcomp,DATAIN,IndicesRB,nDOFsFACEall) ;



TEXTP{end+1} = ['MASTER DOFs = ',num2str(DOFm(:)')] ;

DOFmRB = intersect(IndicesRB,DOFm) ;
TEXTP{end+1} = ['MASTER DOFs (of rigid body type) = ',num2str(DOFmRB(:)')] ;

svdD = svd((Dcomp(DOFm,:)) );
TEXTP{end+1} = ['Ratio SV(end)/SV(1)= ',num2str(svdD(end)/svdD(1)),' (matrix Dcomp_m)'] ;

if rank(Dcomp(DOFm,:)) ~=length(DOFm)
    PrintFileINFO([cd,'/INFO.txt'],TEXTP) ;
    
    error('The selection of Master DOFs is conducive to an ill-conditioned Dcomp(DOFm,:)')
end

if  ~isempty(DOFs)
    TEXTP{end+1} ='WARNING: THIS APPROACH EMPLOYS THE SLAVE DOFS IN ITS FORMULATION';
    [V,DATAOUT] = OutputKinemApproach(Dcomp,DOFs,DOFm,BasisINTFall,nDOFsFACEall,DATAOUT,IndicesRB,BasisINTFcand) ;
    %
else
    % % Number of  DOFs per face
    iacum = 1;
    DOFsFACE = cell(size(nDOFsFACEall)) ;
    DOFsFACE_V = cell(size(nDOFsFACEall)) ;
    iacum_V = 1;
    for iface = 1:length(nDOFsFACEall)
        DOFsFACEloc = iacum:(iacum+nDOFsFACEall(iface)-1) ;
        DOFsFACE{iface} = intersect(DOFsFACEloc,DOFm) ;
        DOFsFACE_V{iface} =  iacum_V:(iacum_V+length(DOFsFACE{iface})-1) ; ;
        iacum = DOFsFACEloc(end)+1 ;
        iacum_V = DOFsFACE_V{iface}(end)+1 ;
    end
    % V.DOFsFACE = DOFsFACE ;
    % V.DOFsFACE_V = DOFsFACE_V ;
    
    % [nDOFsFACE ]= cellfun(@length,DOFsFACE);
    % V.nDOFsFACE = nDOFsFACE ;
    % DATAOUT.BasisInt = V  ;
    
    
    
    V = BasisINTFcand ;
    iacum = 0 ;
    for i = 1:length(V)
        IND_ALL = 1:size(BasisINTFcand{i},2) ;
        IND_ALL = IND_ALL + iacum ;
        [~,INDsel ]= intersect(IND_ALL,DOFsFACE{i}) ;
        V{i} = V{i}(:,INDsel) ;
        V{i} =  full(V{i});
        iacum = iacum + length(IND_ALL) ;
    end
    
    DATAOUT.BasisInt = V  ;
    
end

%
%  p_master = [0 0 0 1 0 0]' ;
% q= Dcomp(DOFm,:)\p_master;
%
% p_slave = Dcomp(DOFs,:)*q ;
%
%  p_slave = Acomp*p_master
%
%  p = zeros(size(Dcomp,1),1) ;
%  p(DOFm ) = p_master ;
%  p(DOFs ) = p_slave ;
%
%  q  =[0 0 0 1 0 0]' ;
%
%  p = Dcomp*q ;
%  q= Dcomp(DOFm,:)\p(DOFm);
%
%  norm(BasisUdom_f*q - BasisINTFcand*p)/norm(BasisUdom_f*q)*100

%
% error('Why is Acomp not full rank ? ')




%
%




% DOFs = 1:size(Dcomp,1) ;
% DOFs = setdiff(DOFs,DOFm) ;
% % In



% for ifgroup=1:length(FACES_GROUPS)
%     f1 = fI{FACES_GROUPS{ifgroup}(1)} ;
%     f2 = fI{FACES_GROUPS{ifgroup}(2)} ;
%     iface = FACES_GROUPS{ifgroup}(1) ;
%     [Vall,RotationMatrixLOC,TEXTP ] = ...
%         ReactionAndInterfaceLocalModes_RVE(BasisUdef,BasisRdef,f1,f2,...
%         DATAIN.TOL_SINGULAR_VALUES_Hqr,...
%         Vrb{iface},M{iface},DATAIN,SinvVal_Udef,SinvVal_Rdef,ifgroup,TEXTP)  ;
%     V{FACES_GROUPS{ifgroup}(1)} = Vall ;
%     V{FACES_GROUPS{ifgroup}(2)} = Vall ;
%     if ~isempty(RotationMatrixLOC)
%         RotationMatrixALL{FACES_GROUPS{ifgroup}(1) } = RotationMatrixLOC{1} ;
%         RotationMatrixALL{FACES_GROUPS{ifgroup}(2) }  = RotationMatrixLOC{2} ;
%     end
% end





%
% V = cell(size(Vrb)) ;
%
% for ifgroup=1:length(FACES_GROUPS)
%     f1 = fI{FACES_GROUPS{ifgroup}(1)} ;
%     f2 = fI{FACES_GROUPS{ifgroup}(2)} ;
%     iface = FACES_GROUPS{ifgroup}(1) ;
%     [Vall,RotationMatrixLOC,TEXTP ] = ...
%         ReactionAndInterfaceLocalModes_RVE(BasisUdef,BasisRdef,f1,f2,...
%         DATAIN.TOL_SINGULAR_VALUES_Hqr,...
%         Vrb{iface},M{iface},DATAIN,SinvVal_Udef,SinvVal_Rdef,ifgroup,TEXTP)  ;
%     V{FACES_GROUPS{ifgroup}(1)} = Vall ;
%     V{FACES_GROUPS{ifgroup}(2)} = Vall ;
%     if ~isempty(RotationMatrixLOC)
%         RotationMatrixALL{FACES_GROUPS{ifgroup}(1) } = RotationMatrixLOC{1} ;
%         RotationMatrixALL{FACES_GROUPS{ifgroup}(2) }  = RotationMatrixLOC{2} ;
%     end
% end
% DATAOUT.BasisInt = V  ;




