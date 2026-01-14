function  [BasisSTRESS,BasisSTRESS_SINGULAR_VALUES ] = ...
    StressModes_viaSVDinel(AsnapSTRESS,DATALOC,DATAoffline,MESH,OPERFE,BasisUdeform,AsnapSTRESSel)  ;

% COMPUTATION OF Stress MODES
% Copy of StressModes_viaSVD.m, explained in 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/104_README_EIFEM_p2D.ml
% In turn, this is the adaptation to the case of inelastic stresses
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
% JAHO, 7-Jan-2024, Balmes, Barcelona 
% ------------------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
  %  DATAoffline.TOLSVD_complementary_modes = 1e-3 ; 
end


% Basic stress modes
% -----------------

% THERE SHOULD BE NO BASIC STRESS MODES HERE 
% disp('-----------------------------------------')
% disp('Basis STRESS modes')
% disp('-----------------------------------------')
%AsnapSTRESS_basic = cell2mat(AsnapSTRESS(DATAoffline.INDEX_BASIC_TESTS)) ;
%     
%  DATALOC.LEGEND_MODES_stress = 'OriginalSnapshotsBasic' ;
% 
% PlotModesStress_SUBDOMAIN(DATALOC,AsnapSTRESS_basic,MESH);


% Basis matrix for
% DATALOC = DefaultField(DATALOC,'TOL',1e-6);
% % -------------------------------------
% 
% % wST = repmat()
% 
% % A = spdiags(Bin,d,6,6);
% 
% % wST = speye(length(OPERFE.wSTs) ; % Diagonal matrix with weights
 nstrain = DATALOC.MESH.nstrain ; 
 wST = repmat(OPERFE.wSTs',nstrain,1);
 wST = wST(:) ; 
 W = spdiags(wST,0,length(wST),length(wST));
% 
% 


%  COMPLEMENTARY MODES STRESSES
disp('-----------------------------------------')
disp('STRESS MODES INELASTIC STRESSES' )
disp('-----------------------------------------')

AsnapSTRESS_compl= cell2mat(AsnapSTRESS(DATAoffline.INDEX_COMPL_TESTS)) ;  % THESE ARE 

if isempty(AsnapSTRESS_compl)  
    disp(['Problem with no inelastic stresses'])
    disp(['Using elastic stresses for computing CECM points (use low tolerance for internal forces if you want to reduce the number of points)'])
   
AsnapSTRESS_compl= cell2mat(AsnapSTRESSel)  ; 

end

% THE STRESSES OBTAINED IN THE COMPLEMENTARY TESTS (PRESUMABLY INELASTIC TESTS)

% 
% AsnapSTRESS_basic= cell2mat(AsnapSTRESSel(DATAoffline.INDEX_COMPL_TESTS)) ; % THESE ARE THE ELASTIC COMPONENT OF THE 
% % ABOVE STRESSES
% 
% 
% [BasisSTRESSbs,S,V,~ ]= WSVDT(AsnapSTRESS_basic,W,DATALOC) ;
% % 
%  disp(['Number of basic stress modes FOR THE INELASTIC TESTS = ',num2str(size(BasisSTRESSbs,2))])
% 
%  DATALOC.LEGEND_MODES_stress = 'Basic' ;
% 
% PlotModesStress_SUBDOMAIN(DATALOC,BasisSTRESSbs,MESH);  





DATALOC = DefaultField(DATALOC,'TOL',1e-6);
 


 
%AsnapSTRESS_compl = AsnapSTRESS_compl/norm(AsnapSTRESS_compl,'fro') ;
%BasisSTRESSbs_ad = BasisSTRESSbs/norm(BasisSTRESSbs,'fro') ;


% Partitioned SVD
Wchol = sqrt(W); 
A = {Wchol*AsnapSTRESS_compl} ;
TOL = DATAoffline.TOLSVD_complementary_modes_stresses; % Tolerance second block
RELTOL = [TOL] ;
DATAlocS = [] ;
[U,S,V] = SRSVD(A,RELTOL,DATAlocS) ;

U = Wchol\U ;  % All modes (including basic ones)

BasisSTRESS = U ; 
BasisSTRESS_SINGULAR_VALUES = S ; 


% 
% = bsxfun(@times,U',S/S(1))' ; 
% BasisSTRESS = bsxfun(@times,U',S/S(1))' ; 


% % ORTHOGONAL COMPLEMENT
% BasisSTRESScomp = SprojDEF_operator(BasisSTRESSbs,W, U) ;
% DATALOCaaa.TOL = 1e-6 ;
% DATALOCaaa.Mchol =Wchol;
% [ BasisSTRESScomp,S,~,~] = WSVDT( BasisSTRESScomp,[],DATALOCaaa) ;

 
 DATALOC.LEGEND_MODES_stress = 'Inelastic_stress' ;

PlotModesStress_SUBDOMAIN(DATALOC,BasisSTRESS,MESH); 

  





% -------------------------------------------------
% POST-PROCESSING SELF-EQUILIBRATED MODES
% -------------------------------------------------
