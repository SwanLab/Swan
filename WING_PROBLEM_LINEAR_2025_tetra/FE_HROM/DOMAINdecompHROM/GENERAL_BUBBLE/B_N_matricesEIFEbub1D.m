function [Bmat,Nmat,WEIGHTSinteg,TRANSF_COORD,CN,DATA,Vrot_all,Kelas,CN_SUPPORT]= B_N_matricesEIFEbub1D(COOR,CN, PROPMAT,MaterialType,DATA,TypeElement,...
    CN_SUPPORT)
% Adaption of B_N_matricesEIFEbub.m for 1D problems
% JAHO, May 13th 2024, Honest Greens, Tuset, Barcelona
% --------------------------------------------------------------
if nargin == 0
    load('tmp2.mat')
end
nelem = size(CN,1);
TRANSF_COORD = cell(nelem,1) ; % INFO ABOUT COORDINATE TRANSFORMATION
Vrot_all = cell(nelem,1) ; % Rotated interface displacements
Bmat = cell(nelem,1) ; % Matrix that maps coarse-scale DOFS onto fine-scale strains at the ECM points of each element
Nmat  = cell(nelem,1) ; % Matrix that maps coarse-scale DOFS onto fine-scale displacements at the ECM points of each element
Kelas = cell(nelem,1) ; %  Elastic stiffness matrix %/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx

WEIGHTSinteg.INTforces = cell(nelem,1) ;  % Determinant of the Jacobian matrix
WEIGHTSinteg.BodyForces = cell(nelem,1) ;  % Determinant of the Jacobian matrix
%DATA =DefaultField(DATA,'IndexPermutationConnectivities',ones(size(MaterialType))) ;
DATA = DefaultField(DATA,'CriterionChooseParentDomain','max_Q_11') ;
DATA = DefaultField(DATA,'UNIFORM_SCALING_REFERENCE_ELEMENT',1) ;   % THIS IS THE DEFAULT OPTION, 11-March-2023
%   Criterion  for choosing the connectivity of the Element  (Cit might be changed within the code in order to minimize the dilatational component
% )


% CECM INTEGRATION RULE FOR NONLINEAR STRESSES (oPTION IMPLEMENTED ON
% 7-Jan-2024), see
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
if isfield(PROPMAT(1).EIFE_prop.INTforces,'CECM_ONLY_FOR_NONLINEAR_STRESSES')
    DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = PROPMAT(1).EIFE_prop.INTforces.CECM_ONLY_FOR_NONLINEAR_STRESSES ;
    if   DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1
    disp('---------- CECM ONLY FOR NONLINEAR PART OF STRESSES -----')
    disp('Coarse-scale stiffness matrix   defined in  PROPMAT(imat).EIFE_prop.Kcoarse  (parent domain ) uses the same elastic constant' )
    disp('as thosed used in training... Do  not change them !  !!!! ')
    end
else
    DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = 0 ;
end

% NOTION OF "EXTENDED" B-MATRICES
% VECTORIZATION IS PREDICATED UPON THE ASSUMPITION THAT ALL NODES HAVE THE
% SAME NUMBER OF DOFS
% THERE, WE HAVE TO DETERMINE WHICH IS THE MAXIMUM NUMBER OF DOFs per node
% Current version (4-oct-2023) assume that the number of CECM points per
% element is the same for all elements
% ------------------------------------------------------------------------

BUBBLE_DOFS = zeros(size(PROPMAT));
for imat = 1:length(PROPMAT)
    BUBBLE_DOFS(imat)= length(PROPMAT(imat).EIFE_prop.INFO.DOFsBUB) ;
    % 26-Dec-2023--->
    
end
nDOFSbub = unique(BUBBLE_DOFS) ;
if length(nDOFSbub) ~= 1
    error('The number of bubble modes is to be the same for all elements in this implementation')
end




for imat = 1:length(PROPMAT)
    PROPMAT(imat) = DefaultField(PROPMAT(imat),'CriterionChooseParentDomain','MAX_Q_11') ;
end

disp('-----------------------------------------------------------')
disp('Determination parent domains + Bmat, Nmat matrices')
disp('-----------------------------------------------------------')


    nnodeE_support = size(CN,2) ;

%DATA.PERMUT = PermutationConnectivities(TypeElement,size(CN,2)) ;
DATA.PERMUT = {[1,2]} ;  


for e = 1:nelem
    disp(['e=',num2str(e)])
    CNloc = CN(e,1:nnodeE_support) ;   %  Nodes of element "e" (ordered by the mesher, in this case GID).
    
    EIFEoper_all = PROPMAT(MaterialType(e)).EIFE_prop ; % Properties EIF object ("element")
  %  DATA.CriterionChooseParentDomain =PROPMAT(MaterialType(e)).CriterionChooseParentDomain ;
    %  [CNnew,TRANSF_COORD{e},Vrot_all{e}] = ParentDomainSearchEIFE(EIFEoper_all,COOR,CNloc,TypeElement,DATA) ;
    [CNnew,TRANSF_COORD{e},Vrot_all{e}] = ParentDomainSearchEIFE1d(EIFEoper_all,COOR,CNloc,TypeElement,DATA) ; % 26-Apr-2024

     CN(e,:) = CNnew ;
   
    % B-MATRIX, weights integration forces
    % --------------------------------------------------------
    [Bmat{e},WEIGHTSinteg.INTforces{e},Kelas{e}] = Bmat_weights_EIFEbub1D(DATA,TRANSF_COORD{e},EIFEoper_all,Vrot_all{e}) ;
    % N-matrix, weights integration body/inertial forces
    [Nmat{e},WEIGHTSinteg.BodyForces{e}] = Nmat_weights_EIFEbub1D(DATA,TRANSF_COORD{e},EIFEoper_all) ;
    
    
    
end

% GHOST DOFS AND EXTENDED VARIABLES
MESHextended = [] ;
[MESHextended.NDOFS_pernode,MESHextended.COOR,MESHextended.CN,MESHextended.DOFS_TO_KEEP,...
    MESHextended.DOFS_bLOC,MESHextended.DOFS_bubLOC,MESHextended.DOFS_ghost] = GhostDOFs(COOR,CN,DATA,nDOFSbub) ;

DATA.MESHextended = MESHextended;

% maxDOF = max(bbb,ndimSP);
% % Maximum number of DOFs per node
% disp(['----------------------------------------------'])
% disp(['Maximum Number of DOFs per node'])
% disp(['DOFmax = ',num2str(maxDOF)])
% disp(['----------------------------------------------'])
% DATA.maxNDOFperNODE = maxDOF ;




