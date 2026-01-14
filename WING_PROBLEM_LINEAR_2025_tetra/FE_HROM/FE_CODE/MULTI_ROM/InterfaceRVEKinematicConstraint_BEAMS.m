function      [V,TEXTP ]= InterfaceRVEKinematicConstraint_BEAMS(Vrb,M,BasisRdef,Vdef_Morth,...
    TOL_SINGULAR_VALUES_Hqr,f1,f2,R,DATAIN,Vrb_Morth,Mdom,BasisUdom,BasisUrb,TEXTP)

%
% [DATAOUT,V,TEXTP] =  InterfaceRVEKinematicConstraint_BEAMS(DATAIN,VrbALL,DATAOUT,...
%     nfaces,fI,BasisUdef,BasisRdef,Mall,...
%     SinvVal_Udef,SinvVal_Rdef,TEXTP,FACES_GROUPS,BasisUrb,Mdom )

if nargin == 0
    load('tmp1.mat')
end
 % Include Kinematic Constraints
% ---------------------------------
nfaces = 2;
BasisINTFdefCAND = cell(nfaces,1) ;
BasisINTFdefCAND{1} = [Vdef_Morth] ;
BasisINTFdefCAND{2} =BasisINTFdefCAND{1} ; % for continuity reasons
VrbALL = cell(nfaces,1) ;
VrbALL{1} = Vrb ;
VrbALL{2} = Vrb ;





% STEP 1. CANDIDATES FOR BEING DEFORMATIONAL INTERFACE MODES  (INTERSECTION)
% ----------------------------------------------
% [BasisINTFdefCAND,BasisRdomROT, BasisUdomROT]= CandidatesInterfaceModes_ANDreact(DATAIN,VrbALL,DATAOUT,...
%     nfaces,fI,BasisUdef,Mall,...
%     TEXTP,FACES_GROUPS,SinvVal_Udef,BasisUrb,Mdom,BasisRdef) ;

% BasisRdom ---> It includes rotation angles




[~,nDEF] = cellfun(@size,BasisINTFdefCAND) ;  %Number of deformational modes
[~,nRB] = cellfun(@size,VrbALL) ;   % Number of RB modes

nCAND = nDEF(:) + nRB(:) ;

% STEP 2. All candidates
% ----------------------
BasisINTFcand =  cell(size(BasisINTFdefCAND)) ;
BasisINTFcandROT =  cell(size(BasisINTFdefCAND)) ;


 

IndicesRB = [] ;
%IndicesDEF = [] ;
iacum = 0 ;
for iface = 1:length(VrbALL)
    IndicesRBloc =  (1:nRB(iface))+iacum ;
    
    IndicesRB = [IndicesRB,IndicesRBloc] ;
    BasisINTFcand{iface} = sparse([VrbALL{iface},BasisINTFdefCAND{iface}]) ;
    
  
    
    %     IndicesDEFloc =  (nRB(iface):size(BasisINTFcand{iface},2)  )+iacum ;
    %      IndicesDEF{end+1} =
    iacum = iacum +  nRB(iface) + nDEF(iface) ;
end
DATAIN = DefaultField(DATAIN,'ROTATION_LOC',[]) ;
R = DATAIN.ROTATION_LOC ;
BasisINTFcandROT=  BasisINTFcand ;
if ~isempty(R)
    BasisINTFcandROT{2} =  RotateMatrix(R,    BasisINTFcandROT{2} )  ;
end





[nelems,nDOFsFACEall] = cellfun(@size,BasisINTFcand) ;

IndicesDEF = 1:sum(nCAND) ;
IndicesDEF(IndicesRB) = [] ;

% STEP 3. Form a diagonal basis matrix
% ------------------------------------
BasisINTFall = blkdiag(BasisINTFcand{:}) ;
BasisINTFallROT = blkdiag(BasisINTFcandROT{:}) ;

% STEP 4. Basis of domain displacements
% f = cell2mat(fI') ;
% BasisUrb_f =BasisUrb(f,:) ;
% BasisUdef_f = BasisUdef(f,:) ;
% % % RB modes with norm = 1
% for i = 1:size(BasisUrb_f,2)
%     BasisUrb_f(:,i) = BasisUrb_f(:,i)/norm(BasisUrb_f(:,i)) ;
% end
% for i = 1:size(BasisUdef_f,2)
%     BasisUdef_f(:,i) = BasisUdef_f(:,i)/norm(BasisUdef_f(:,i))*SinvVal_Udef(i)/SinvVal_Udef(1) ;
% end

%BasisUdom_f = [BasisUrb_f,BasisUdef_f] ;

% STEP 5. Solving minimization problem
% 
% 
 DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,...
     'WORK_DONE_BY_REACTIVE_FORCES',0);
    f = [f1;f2] ; 

 if DATAIN.KINEMATIC_CONSTRAINTS_MODES.WORK_DONE_BY_REACTIVE_FORCES == 1
%     % May-14th --- Approach based on the work done by reactive forces
BasisRdom = [Mdom*BasisUrb,BasisRdef] ; 
     Dcomp = BasisINTFallROT'*BasisRdom(f,:);
%     
 else
    % This is equivalent to solve a least-squares problem in the norm of
    % the mass matrix. Pertinent tests were made and the results appeared
    % to confirm that the implementation is in principle correct
    Mchol =   chol(Mdom(f,f)) ;
    MBasisINTFall = Mchol*BasisINTFallROT ;
    % MBasisUdom = Mchol*BasisUdom_f ;
 %   BasisUdom = cell2mat(BasisUdom) ; 
    MBasisUdom = Mchol*BasisUdom(f,:) ;
    
    
    Dcomp = MBasisINTFall\MBasisUdom ;
end

%    Dcomp = BasisINTFcand\BasisUdom_f ;

% if mod(size(BasisUdom_f,2),2)~=0
%     error('The total number of domain modes must be even')
% end

% STEP 6. Master and Slave DOFs
% Select linearly independent columns
% --------------------------------------
% DATAIN.KINEMATIC_CONSTRAINTS_MODES = ...
%     DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'DEIM_BASED_METHOD',1) ;

%if DATAIN.KINEMATIC_CONSTRAINTS_MODES.DEIM_BASED_METHOD == 1
    [DOFm,DOFs] = CoarseDOFS_DEIMbased(Dcomp,DATAIN,IndicesRB,nDOFsFACEall) ;
    
    
%else
    % Before 10-May-2019
 %   [DOFm,DOFs] = CoarseDOFS_master_slave(Dcomp,DATAIN,BasisUdom_f,IndicesRB,IndicesDEF,nDOFsFACEall,nRB) ;
    
%end

TEXTP{end+1} = ['MASTER DOFs = ',num2str(DOFm(:)')] ;

DOFmRB = intersect(IndicesRB,DOFm) ;
TEXTP{end+1} = ['MASTER DOFs (of rigid body type) = ',num2str(DOFmRB(:)')] ;

svdD = svd((Dcomp(DOFm,:)) );
TEXTP{end+1} = ['Ratio SV(end)/SV(1)= ',num2str(svdD(end)/svdD(1)),' (matrix Dcomp_m)'] ;

if rank(Dcomp(DOFm,:)) ~=length(DOFm)
    PrintFileINFO([cd,'/INFO.txt'],TEXTP) ;
    
    error('The selection of Master DOFs is conducive to an ill-conditioned Dcomp(DOFm,:)')
end

%Acomp = Dcomp(DOFs,:)*inv(Dcomp(DOFm,:)) ;

% OUTPUT WILL BE A STRUCTURED ARRAY

%V.DOFmP = DOFm ;   % Master DOFs
%V.DOFsP = DOFs ;   % Slave DOFs
%V.IndicesRB = IndicesRB ; % Indexes Rigid Body DOFs
%V.BasisINTFall_cell = BasisINTFcand ;  % Basis containing all modes (cell array)
BasisINT = BasisINTFall(:,DOFm)  ; % + BasisINTFall(:,DOFs)*Acomp ;   % Basis matrix to be used in the formulatiion
nDOF_1 = length(DOFm)/2 ; 

nDOF_rows = size(BasisINT,1)/2 ; 

V  = full(BasisINT(1:nDOF_rows,1:nDOF_1)) ; 
 
 

%V.Acomp = Acomp ;
% Number of  DOFs per face
% iacum = 1;
% DOFsFACE = cell(size(nDOFsFACEall)) ;
% DOFsFACE_V = cell(size(nDOFsFACEall)) ;
% iacum_V = 1;
% for iface = 1:length(nDOFsFACEall)
%     DOFsFACEloc = iacum:(iacum+nDOFsFACEall(iface)-1) ;
%     DOFsFACE{iface} = intersect(DOFsFACEloc,DOFm) ;
%     DOFsFACE_V{iface} =  iacum_V:(iacum_V+length(DOFsFACE{iface})-1) ; ;
%     iacum = DOFsFACEloc(end)+1 ;
%     iacum_V = DOFsFACE_V{iface}(end)+1 ;
% end
% V.DOFsFACE = DOFsFACE ;
% V.DOFsFACE_V = DOFsFACE_V ;
% 
% [nDOFsFACE ]= cellfun(@length,DOFsFACE);
% V.nDOFsFACE = nDOFsFACE ;
% DATAOUT.BasisInt = V  ;



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
