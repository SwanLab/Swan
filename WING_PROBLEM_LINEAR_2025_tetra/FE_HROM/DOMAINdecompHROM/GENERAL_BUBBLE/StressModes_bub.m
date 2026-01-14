function  [BasisSTRESS,BasisSTRESS_SINGULAR_VALUES ] =   StressModes_plas(AsnapSTRESS,DATALOC,DATAoffline,MESH,OPERFE,BasisUdeform)  ;

% COMPUTATION OF Stress MODES
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/104_README_EIFEM_p2D.mlx


if nargin == 0
    load('tmp1.mat')
    
end


% % Basic stress modes
% % -----------------
% disp('-----------------------------------------')
% disp('Basis STRESS modes')
% disp('-----------------------------------------')
% % % 
% % % % -----------------------------------------------------------------------
%   SNAPbasic = cell2mat(AsnapSTRESS(DATAoffline.INDEX_BASIC_TESTS)) ;
%  SNAPcompl = AsnapSTRESS(DATAoffline.INDEX_COMPL_TESTS) ;



%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% disp('----------------------------------------------------------------------')
 disp(['Computing  modes via right-singular vectors displacements'])
% disp('----------------------------------------------------------------------')
% % A = U*S*V^T -->   U = (A*V)*inv(S)
% 
% disp(['This operation may be optimized...'])
AsnapSTRESS = cell2mat(AsnapSTRESS) ; 
 BasisSTRESS = AsnapSTRESS*BasisUdeform.V_upsilonDEF;
% 
%   

% wST = speye(length(OPERFE.wSTs) ; % Diagonal matrix with weights
nstrain = DATALOC.MESH.nstrain ; 
wST = repmat(OPERFE.wSTs',nstrain,1);
wST = wST(:) ; 
W = spdiags(wST,0,length(wST),length(wST));

DATALOCaaa = [] ; 
[BasisSTRESS,S,V,~ ]= WSVDT(BasisSTRESS,W,DATALOCaaa) ;

BasisSTRESS_SINGULAR_VALUES = S/S(1); 

disp(['Number of   stress modes = ',num2str(size(BasisSTRESS,2))])

 DATALOC.LEGEND_MODES_stress = 'STRESS' ;
DATALOC.NAME_BASE = '' ; 
PlotModesStress_SUBDOMAIN(DATALOC,BasisSTRESS,MESH);  

CHECK_error = 1; 

if CHECK_error == 1
   BasisSTRESSC = SVDT(BasisSTRESS) ; 
  
    
    c = BasisSTRESSC'*AsnapSTRESS ;
    EEE = norm(AsnapSTRESS-BasisSTRESSC*c,'fro') ; 
    nn = norm(AsnapSTRESS ,'fro') ;
    disp(['ERROR_STRESSES  (over 1)=',num2str(EEE/nn)])
    
     
end




% 
% 
% 
% %  COMPLEMENTARY MODES STRESSES
% disp('-----------------------------------------')
% disp('ALL STRESS modes' )
% disp('-----------------------------------------')
% 
% AsnapSTRESS_compl= cell2mat(AsnapSTRESS(DATAoffline.INDEX_COMPL_TESTS)) ;
%  
% AsnapSTRESS_compl = AsnapSTRESS_compl/norm(AsnapSTRESS_compl,'fro') ;
% BasisSTRESSbs_ad = BasisSTRESSbs/norm(BasisSTRESSbs,'fro') ;
% 
% 
% % Partitioned SVD
% Wchol = sqrt(W); 
% A = {Wchol*BasisSTRESSbs_ad,Wchol*AsnapSTRESS_compl} ;
% TOL = DATAoffline.TOLSVD_complementary_modes_stresses; % Tolerance second block
% RELTOL = [0,TOL] ;
% DATAlocS = [] ;
% [U,S,V] = SRSVD(A,RELTOL,DATAlocS) ;
% 
% U = Wchol\U ;  % All modes (including basic ones)
% 
% BasisSTRESS = U ; 
% BasisSTRESS_SINGULAR_VALUES = S ; 
% 
% 
% % 
% % = bsxfun(@times,U',S/S(1))' ; 
% % BasisSTRESS = bsxfun(@times,U',S/S(1))' ; 
% 
% 
% % % ORTHOGONAL COMPLEMENT
% % BasisSTRESScomp = SprojDEF_operator(BasisSTRESSbs,W, U) ;
% % DATALOCaaa.TOL = 1e-6 ;
% % DATALOCaaa.Mchol =Wchol;
% % [ BasisSTRESScomp,S,~,~] = WSVDT( BasisSTRESScomp,[],DATALOCaaa) ;
% 
%  
%  DATALOC.LEGEND_MODES_stress = 'All_stresses' ;
% 
% PlotModesStress_SUBDOMAIN(DATALOC,BasisSTRESS,MESH); 

%   
% 
% 
% 
% 
% 
% % -------------------------------------------------
% % POST-PROCESSING SELF-EQUILIBRATED MODES
% % -------------------------------------------------
