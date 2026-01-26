function [Bmat,Nmat,WEIGHTSinteg,TRANSF_COORD,CN,DATA,Vrot_all]= B_N_matricesEIFE(COOR,CN, PROPMAT,MaterialType,DATA,TypeElement) 
% EIFE METHOD, DETERMINATION OF B and N matrices  --------------
% --------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;
TRANSF_COORD = cell(nelem,1) ; % INFO ABOUT COORDINATE TRANSFORMATION
Vrot_all = cell(nelem,1) ; % Rotated interface displacements
Bmat = cell(nelem,1) ; % Matrix that maps coarse-scale DOFS onto fine-scale strains at the ECM points of each element
Nmat  = cell(nelem,1) ; % Matrix that maps coarse-scale DOFS onto fine-scale displacements at the ECM points of each element
WEIGHTSinteg.INTforces = cell(nelem,1) ;  % Determinant of the Jacobian matrix
WEIGHTSinteg.BodyForces = cell(nelem,1) ;  % Determinant of the Jacobian matrix
%DATA =DefaultField(DATA,'IndexPermutationConnectivities',ones(size(MaterialType))) ;
DATA = DefaultField(DATA,'CriterionChooseParentDomain','max_Q_11') ;
DATA = DefaultField(DATA,'UNIFORM_SCALING_REFERENCE_ELEMENT',1) ;   % THIS IS THE DEFAULT OPTION, 11-March-2023

%   Criterion  for choosing the connectivity of the Element  (Cit might be changed within the code in order to minimize the dilatational component
% )
for imat = 1:length(PROPMAT)
    PROPMAT(imat) = DefaultField(PROPMAT(imat),'CriterionChooseParentDomain','MAX_Q_11') ;
end

disp('-----------------------------------------------------------')
disp('Determination parent domains + Bmat, Nmat matrices')
disp('-----------------------------------------------------------')
DATA.PERMUT = PermutationConnectivities(TypeElement,size(CN,2)) ; 
for e = 1:nelem    
    disp(['e=',num2str(e)])
    CNloc = CN(e,:) ;   %  Nodes of element "e" (ordered by the mesher, in this case GID)
    EIFEoper_all = PROPMAT(MaterialType(e)).EIFE_prop ; % Properties EIF object ("element")
    DATA.CriterionChooseParentDomain =PROPMAT(MaterialType(e)).CriterionChooseParentDomain ;
    [CNnew,TRANSF_COORD{e},Vrot_all{e}] = ParentDomainSearchEIFE(EIFEoper_all,COOR,CNloc,TypeElement,DATA) ;
    CN(e,:) = CNnew ;   % Coordinates of the nodes of element "e"    
    % B-MATRIX, weights integration forces
    % --------------------------------------------------------
    [Bmat{e},WEIGHTSinteg.INTforces{e}] = Bmat_weights_EIFE(DATA,TRANSF_COORD{e},EIFEoper_all,Vrot_all{e}) ;     
    % N-matrix, weights integration body/inertial forces 
    [Nmat{e},WEIGHTSinteg.BodyForces{e}] = Nmat_weights_EIFE(DATA,TRANSF_COORD{e},EIFEoper_all) ; 
    
    
    
end




