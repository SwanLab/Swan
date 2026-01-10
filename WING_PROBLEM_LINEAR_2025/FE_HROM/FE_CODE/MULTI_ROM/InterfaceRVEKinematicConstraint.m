function  [DATAOUT,V,TEXTP] =  InterfaceRVEKinematicConstraint(DATAIN,VrbALL,DATAOUT,...
    nfaces,fI,BasisUdef,BasisRdef,Mall,...
    SinvVal_Udef,SinvVal_Rdef,TEXTP,FACES_GROUPS,BasisUrb,Mdom )

if nargin == 0
    load('tmp1.mat')
   % DATAIN.KINEMATIC_CONSTRAINTS_MODES.TOL_INTERSECT_DISP = 1;   %  Tolerance ANGLE (deg.) determining intersection pair of faces
   % DATAIN.KINEMATIC_CONSTRAINTS_MODES.TOL_SVD_DISP_MODES = 1e-3;   %  Tolerance for the SVD used in computing the intersection
    
end
% Include Kinematic Constraints
% ---------------------------------

 
% STEP 1. CANDIDATES FOR BEING DEFORMATIONAL INTERFACE MODES  (INTERSECTION)
% ----------------------------------------------
% [BasisINTFdefCAND,BasisRdomROT, BasisUdomROT]= CandidatesInterfaceModesW_RVE(DATAIN,VrbALL,DATAOUT,...
%     nfaces,fI,BasisUdef,Mall,...
%     TEXTP,FACES_GROUPS,SinvVal_Udef,BasisUrb,DATA_REFMESH.M,BasisRdef) ;

DATAIN.KINEMATIC_CONSTRAINTS_MODES = ...
    DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'FILTERING_DISP_USING_WORK_REACTIONS',0) ;

if DATAIN.KINEMATIC_CONSTRAINTS_MODES.FILTERING_DISP_USING_WORK_REACTIONS  == 0
    
    [BasisINTFdefCAND,BasisRdomROT, BasisUdomROT]= CandidatesInterfaceModes_ANDreact(DATAIN,VrbALL,DATAOUT,...
        nfaces,fI,BasisUdef,Mall,...
        TEXTP,FACES_GROUPS,SinvVal_Udef,BasisUrb,Mdom,BasisRdef) ;
    BasisINTFrbORTH = VrbALL ; 
else
    
    [BasisINTFdefCAND,BasisRdomROT, BasisUdomROT,TEXTP,BasisINTFrbORTH]= ...
        CandidatesInterfaceModes_WCOMB(DATAIN,VrbALL,DATAOUT,...
        nfaces,fI,BasisUdef,Mall,...
        TEXTP,FACES_GROUPS,SinvVal_Udef,BasisUrb,Mdom,BasisRdef,SinvVal_Rdef) ;
    
end

% BasisRdom ---> It includes rotation angles 




[~,nDEF] = cellfun(@size,BasisINTFdefCAND) ;  %Number of deformational modes
[~,nRB] = cellfun(@size,VrbALL) ;   % Number of RB modes

nCAND = nDEF(:) + nRB(:) ;

% STEP 2. All candidates
% ----------------------
BasisINTFcand =  cell(size(BasisINTFdefCAND)) ;
BasisINTFcandORTH = cell(size(BasisINTFdefCAND)) ; 
IndicesRB = [] ;
%IndicesDEF = [] ;
iacum = 0 ;
for iface = 1:length(VrbALL)
    IndicesRBloc =  (1:nRB(iface))+iacum ;
    
    IndicesRB = [IndicesRB,IndicesRBloc] ;
    BasisINTFcand{iface} = sparse([VrbALL{iface},BasisINTFdefCAND{iface}]) ;
    BasisINTFcandORTH{iface} = [BasisINTFrbORTH{iface},BasisINTFdefCAND{iface}]
    %     IndicesDEFloc =  (nRB(iface):size(BasisINTFcand{iface},2)  )+iacum ;
    %      IndicesDEF{end+1} =
    iacum = iacum +  nRB(iface) + nDEF(iface) ;
end
[nelems,nDOFsFACEall] = cellfun(@size,BasisINTFcand) ;

IndicesDEF = 1:sum(nCAND) ;
IndicesDEF(IndicesRB) = [] ;

% STEP 3. Form a diagonal basis matrix
% ------------------------------------
BasisINTFall = blkdiag(BasisINTFcand{:}) ;

% STEP 4. Basis of domain displacements
f = cell2mat(fI') ;
BasisUrb_f =BasisUrb(f,:) ;
BasisUdef_f = BasisUdef(f,:) ;
% % % RB modes with norm = 1
% for i = 1:size(BasisUrb_f,2)
%     BasisUrb_f(:,i) = BasisUrb_f(:,i)/norm(BasisUrb_f(:,i)) ;
% end
% for i = 1:size(BasisUdef_f,2)
%     BasisUdef_f(:,i) = BasisUdef_f(:,i)/norm(BasisUdef_f(:,i))*SinvVal_Udef(i)/SinvVal_Udef(1) ;
% end

BasisUdom_f = [BasisUrb_f,BasisUdef_f] ;

% Temporary. Criterion to filter out candidate modes
% Rotations are ignored. May-19-2019
% T{1} = BasisRdef(fI{1},:)'*BasisINTFcandORTH{1} ; % + BasisRdef(fI{3},:)'*BasisINTFcandORTH{3} ; 
% T{2} = BasisRdef(fI{2},:)'*BasisINTFcandORTH{2} ; % + BasisRdef(fI{4},:)'*BasisINTFcandORTH{4} ; 
% T{3} = BasisRdef(fI{3},:)'*BasisINTFcandORTH{3} ; % + BasisRdef(fI{4},:)'*BasisINTFcandORTH{4} ; 
% T{4} = BasisRdef(fI{4},:)'*BasisINTFcandORTH{4} ; 
% 
% S12 = svd([T{1} T{2}]) ; 
% S34 = svd([T{3} T{4}]) ; 
% S123 = svd([T{1} T{2} T{3}]) ; 
% S13 = svd([T{1} T{3}]) 

BasisUDEFf12 =[  BasisUdef([fI{1}; fI{2}],:)];
 
% -----------------------------------------------



% STEP 5. Solving minimization problem
if size(BasisINTFall,2) <= size(BasisUdom_f)
    error('Disable option of kinematical constraints')
else
    
 
    DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,...
        'WORK_DONE_BY_REACTIVE_FORCES',0); 

    if DATAIN.KINEMATIC_CONSTRAINTS_MODES.WORK_DONE_BY_REACTIVE_FORCES == 1
        % May-14th --- Approach based on the work done by reactive forces 
         Dcomp = BasisINTFall'*BasisRdomROT ; 
        
    else
    % This is equivalent to solve a least-squares problem in the norm of
    % the mass matrix. Pertinent tests were made and the results appeared
    % to confirm that the implementation is in principle correct
  
    Mchol =   chol(Mdom(f,f)) ;
    MBasisINTFall = Mchol*BasisINTFall ;
   % MBasisUdom = Mchol*BasisUdom_f ;
        MBasisUdom = Mchol*BasisUdomROT ;

    
    Dcomp = MBasisINTFall\MBasisUdom ;
    end
    
    %    Dcomp = BasisINTFcand\BasisUdom_f ;
end

% if mod(size(BasisUdom_f,2),2)~=0
%     error('The total number of domain modes must be even')
% end

% STEP 6. Master and Slave DOFs
% Select linearly independent columns
% --------------------------------------
DATAIN.KINEMATIC_CONSTRAINTS_MODES = ...
    DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'DEIM_BASED_METHOD',1) ; 

if DATAIN.KINEMATIC_CONSTRAINTS_MODES.DEIM_BASED_METHOD == 1
    [DOFm,DOFs] = CoarseDOFS_DEIMbased(Dcomp,DATAIN,IndicesRB,nDOFsFACEall) ;

    
else
    % Before 10-May-2019
[DOFm,DOFs] = CoarseDOFS_master_slave(Dcomp,DATAIN,BasisUdom_f,IndicesRB,IndicesDEF,nDOFsFACEall,nRB) ;

end

TEXTP{end+1} = ['MASTER DOFs = ',num2str(DOFm(:)')] ;

DOFmRB = intersect(IndicesRB,DOFm) ;
TEXTP{end+1} = ['MASTER DOFs (of rigid body type) = ',num2str(DOFmRB(:)')] ;

svdD = svd((Dcomp(DOFm,:)) );
TEXTP{end+1} = ['Ratio SV(end)/SV(1)= ',num2str(svdD(end)/svdD(1)),' (matrix Dcomp_m)'] ;

if rank(Dcomp(DOFm,:)) ~=length(DOFm)
    PrintFileINFO([cd,'/INFO.txt'],TEXTP) ;
    
    error('The selection of Master DOFs is conducive to an ill-conditioned Dcomp(DOFm,:)')
end


% CASE YOU INCLUDE SLAVE DOFS
[V,DATAOUT] = OutputKinemApproach(Dcomp,DOFs,DOFm,BasisINTFall,nDOFsFACEall,DATAOUT,IndicesRB,BasisINTFcand) ; 
 
% CASE YOU DON'T INCLUDE SLAVE DOFS

if isempty(DOFs)
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
