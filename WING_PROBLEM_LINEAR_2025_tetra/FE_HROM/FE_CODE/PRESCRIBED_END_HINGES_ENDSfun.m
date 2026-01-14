function  [G,uBAR,DOFr,DOFm,AREA,R,DOFA,DOFB] = ...
    PRESCRIBED_END_HINGES_ENDSfun(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)

% See BeamROM_nonlinear.pdf, /home/joaquin/Desktop/CURRENT_TASKS/POTENTIAL_RESEARCH_TOPICS/COMBINING_MULTISCALE_REDUCTIONMODELS/REPORT_MULTIS_REDUC_MODEL/PAPER_multiLEARN 
% General boundary conditions --- Free warping 

if nargin == 0 
    load('tmp1.mat')
end

if ~iscell(a_A)
    a_A = a_A(:) ; a_A = a_A' ; 
    a_A = mat2cell(a_A,1,ones(1,length(a_A))) ;   
end
if ~iscell(a_B)
    a_B = a_B(:) ; a_B = a_B' ; 
    a_B = mat2cell(a_B,1,ones(1,length(a_B))) ;   
end
% Restrictions each end 
[s_A] = cellfun(@isempty,a_A) ;
[s_B] = cellfun(@isempty,a_B) ;
s_A = find(s_A == 0) ; 
s_B = find(s_B == 0) ; 
a_A = cell2mat(a_A(s_A)) ;
a_B = cell2mat(a_B(s_B)) ;
a_As = a_A(:) ; 
a_Bs = a_B(:) ; 

nmodes = 6; 
if size(COOR,2) ==2
    nmodes= 3; 
end

 




ndim = size(COOR,2); 
%%%% FACE 1
iface=1 ; 
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ; 
if  max(nodesfA-sort(nodesfA)) ~=0 
    error(['nodesfA should be sorted in ascending order'])
end
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
% CentoidFA = sum(COOR_FACE,1)/size(COOR_FACE,1); % Centroid
% COORrelA = bsxfun(@minus,COOR_FACE',CentoidFA')'; % Coordinates relative to centroid
% R = ConstructBasisRigidBody(COORrelA) ; % Basis matrix for rigid body motions, relative to centroid
%%% Renumbering 
% Local coordinates  --> COOR_FACE 
% Local connectivities 
[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{iface},TypeElementB) ; 


% [Nst,wST] = GeometricMassMatrixSurface(CONNECTb{1},nodesfA,COOR_FACE,TypeElementB) ; 
% wSTdiag = CompWeightDiag(wST,1)  ; 
% Mst = (wSTdiag*Nst)'*Nst ; 
% % Recomputing centroid 
% CentroidFA = zeros(1,3) ; 
% AREA = sum(wST) ; 
% for idim = 1:3 
%     CentroidFA(idim) = wST'*(Nst*COOR_FACE(:,idim))/AREA ; 
% end
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ;

 
% FACE 2
iface=2 ; 
%nodesf = unique(CONNECTb{end,iface}) ;    % Nodes face2
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired, 
 
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B

% Rorth = R'*M --->  R(INDX)'*N'*W*N 

Rorth = zeros(size(R')) ; 
for idim =1:ndim
    INDLOC =idim:ndim:size(R,1) ; 
    Rorth(:,INDLOC) = R(INDLOC,:)'*Mst ; 
end
Rorth = Rorth'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
 
% -------- 
[JA,bA,rA,lA] = CondFaceLoc(R,s_A,Rorth,a_As,DOFA) ; 
[JB,bB,rB,lB] = CondFaceLoc(R,s_B,Rorth,a_Bs,DOFB) ; 
 
% 
DOFr = [DOFA(rA); DOFB(rB) ]  ; % Master DOFs
DOFm = [DOFA(lA); DOFB(lB) ]  ; % Slave DOFs 
% ----------------------------

G = sparse(length(DOFr),length(DOFm)) ; 
uBAR = zeros(length(DOFr),1) ; 

% -A 
iini = 1; 
ifin = size(JA,1) ; 
jini = 1; 
jfin = size(JA,2) ; 

G(iini:ifin,jini:jfin) = JA ; 
uBAR(iini:ifin) = bA  ; 
% -B
iini = ifin+1; 
ifin = length(DOFr) ; 
jini = jfin+1; 
jfin = length(DOFm) ; 

G(iini:ifin,jini:jfin) = JB; 
uBAR(iini:ifin) = bB  ; 





end 


function  [JA,bA,rA,lA] = CondFaceLoc(R,s_A,Rorth,a_As,DOFA)

R_As = R(:,s_A) ; 
Rbar_As = Rorth(:,s_A) ; 
% Face A 
[~,rA]=licols(Rbar_As') ; %
lA = setdiff(1:length(DOFA),rA) ;
JA = -Rbar_As(rA,:)'\Rbar_As(lA,:)' ;
bA = Rbar_As(rA,:)'\(Rbar_As'*R_As*a_As   ) ; 

end

 

 