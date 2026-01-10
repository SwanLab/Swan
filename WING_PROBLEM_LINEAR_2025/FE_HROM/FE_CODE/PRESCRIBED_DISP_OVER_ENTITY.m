function  [G,uBAR,DOFr,DOFm ] = ...
    PRESCRIBED_DISP_OVER_ENTITY(a_A,nodesfA,COOR,CONNECTb,TypeElementB,...
    DATA)

% See BeamROM_nonlinear.pdf, /home/joaquin/Desktop/CURRENT_TASKS/POTENTIAL_RESEARCH_TOPICS/COMBINING_MULTISCALE_REDUCTIONMODELS/REPORT_MULTIS_REDUC_MODEL/PAPER_multiLEARN 
% General boundary conditions ---  Copy of PRESCRIBED_END_HINGES_ENDSfun.m 

if nargin == 0 
    load('tmp2.mat')
end



 
% Restrictions each end
[s_A] = cellfun(@isempty,a_A) ;

 
s_A = find(s_A == 0) ;
a_A = cell2mat(a_A(s_A)) ;
a_As = a_A(:) ;

nmodes = 6; ndim = 3; 
if size(COOR,2) ==2
    nmodes= 3; 
    ndim = 2; 
end

 
 
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
% CentoidFA = sum(COOR_FACE,1)/size(COOR_FACE,1); % Centroid
% COORrelA = bsxfun(@minus,COOR_FACE',CentoidFA')'; % Coordinates relative to centroid
% R = ConstructBasisRigidBody(COORrelA) ; % Basis matrix for rigid body motions, relative to centroid
%%% Renumbering 
% Local coordinates  --> COOR_FACE 
% Local connectivities 

[CNb setBelem]= ElemBnd(CONNECTb,nodesfA) ; 

[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CNb,TypeElementB) ; 

 
% end
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ;

  

% Rorth = R'*M --->  R(INDX)'*N'*W*N 

Rorth = zeros(size(R')) ; 
for idim =1:ndim
    INDLOC =idim:ndim:size(R,1) ; 
    Rorth(:,INDLOC) = R(INDLOC,:)'*Mst ; 
end
Rorth = Rorth'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
if length(s_A) == nmodes 
    % All  DOFs  are prescribed 
    DOFr = DOFA; 
    DOFm = [] ; 
    G = [] ; 
    uBAR = R*a_A; 
    
    
else
 
[JA,bA,rA,lA] = CondFaceLoc(R,s_A,Rorth,a_As,DOFA) ;  
% 
DOFr = [DOFA(rA)  ]  ; % Master DOFs
DOFm = [DOFA(lA)]  ; % Slave DOFs 
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
% % -B
% iini = ifin+1; 
% ifin = length(DOFr) ; 
% jini = jfin+1; 
% jfin = length(DOFm) ; 
% G(iini:ifin,jini:jfin) = JB; 
% uBAR(iini:ifin) = bB  ; 

end 

end


function  [JA,bA,rA,lA] = CondFaceLoc(R,s_A,Rorth,a_As,DOFA)

R_As = R(:,s_A) ; 
Rbar_As = Rorth(:,s_A) ; 
% Face A 
[~,rA]=licols(Rbar_As') ; %
lA = setdiff(1:length(DOFA),rA) ;
JA = -Rbar_As(rA,:)'\Rbar_As(lA,:)' ;
bA = Rbar_As(rA,:)'\(Rbar_As'*R_As*a_As   ) ; 
JA = sparse(JA) ; 

end

 

 