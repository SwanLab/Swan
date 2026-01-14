function  [BasisSTRESS,BasisSTRESS_SINGULAR_VALUES ] =   StressModes_viaSVD(AsnapSTRESS,DATALOC,DATAoffline,MESH,OPERFE,BasisUdeform)  ;

% COMPUTATION OF Stress MODES
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/104_README_EIFEM_p2D.ml
if nargin == 0
    load('tmp.mat')
    DATAoffline.TOLSVD_complementary_modes = 1e-3 ; 
end


% Basic stress modes
% -----------------
disp('-----------------------------------------')
disp('Basis STRESS modes')
disp('-----------------------------------------')

 


AsnapSTRESS_basic = cell2mat(AsnapSTRESS(DATAoffline.INDEX_BASIC_TESTS)) ;

%     
%  DATALOC.LEGEND_MODES_stress = 'OriginalSnapshotsBasic' ;
% 
% PlotModesStress_SUBDOMAIN(DATALOC,AsnapSTRESS_basic,MESH);


% Basis matrix for
DATALOC = DefaultField(DATALOC,'TOL',1e-6);
% -------------------------------------

% wST = repmat()

% A = spdiags(Bin,d,6,6);

% wST = speye(length(OPERFE.wSTs) ; % Diagonal matrix with weights
nstrain = DATALOC.MESH.nstrain ; 
ISLARGE = 0 ; 
if size(AsnapSTRESS_basic,1) ~= DATALOC.MESH.ndofSTRESS
    % Large strain formulation 
    nstrain = DATALOC.MESH.ndim^2; 
    ISLARGE = 1; 
end

wST = repmat(OPERFE.wSTs',nstrain,1);
wST = wST(:) ; 
W = spdiags(wST,0,length(wST),length(wST));


[BasisSTRESSbs,S,V,~ ]= WSVDT(AsnapSTRESS_basic,W,DATALOC) ;

disp(['Number of basic stress modes = ',num2str(size(BasisSTRESSbs,2))])


PLOT_BASIS_MODES = 0 ; 

if PLOT_BASIS_MODES == 1
 DATALOC.LEGEND_MODES_stress = 'Basic' ;

PlotModesStress_SUBDOMAIN(DATALOC,BasisSTRESSbs,MESH);  

end

%  COMPLEMENTARY MODES STRESSES
disp('-----------------------------------------')
disp('ALL STRESS modes' )
disp('-----------------------------------------')

AsnapSTRESS_compl= cell2mat(AsnapSTRESS(DATAoffline.INDEX_COMPL_TESTS)) ;
 
AsnapSTRESS_compl = AsnapSTRESS_compl/norm(AsnapSTRESS_compl,'fro') ;
BasisSTRESSbs_ad = BasisSTRESSbs/norm(BasisSTRESSbs,'fro') ;


% Partitioned SVD
Wchol = sqrt(W); 
if ~isempty(AsnapSTRESS_compl)
A = {Wchol*BasisSTRESSbs_ad,Wchol*AsnapSTRESS_compl} ;
TOL = DATAoffline.TOLSVD_complementary_modes_stresses; % Tolerance second block
RELTOL = [0,TOL] ;
DATAlocS = [] ;

else
A = {Wchol*BasisSTRESSbs_ad} ;
RELTOL = [0] ;
DATAlocS = [] ;
end
[U,S,V] = SRSVD(A,RELTOL,DATAlocS) ;

U = Wchol\U ;  % All modes (including basic ones)

BasisSTRESS = U ; 
BasisSTRESS_SINGULAR_VALUES = S/S(1) ; 


% 
% = bsxfun(@times,U',S/S(1))' ; 
% BasisSTRESS = bsxfun(@times,U',S/S(1))' ; 


% % ORTHOGONAL COMPLEMENT
% BasisSTRESScomp = SprojDEF_operator(BasisSTRESSbs,W, U) ;
% DATALOCaaa.TOL = 1e-6 ;
% DATALOCaaa.Mchol =Wchol;
% [ BasisSTRESScomp,S,~,~] = WSVDT( BasisSTRESScomp,[],DATALOCaaa) ;

 PLOT_BASIS_MODES = 0 ; 

if PLOT_BASIS_MODES == 1

 disp('PLotting stress modes... ')
if ISLARGE == 1
 disp('Warning = PK1 stresses are not correctly post-process in GID (they are not symmetric)')
end

 DATALOC.LEGEND_MODES_stress = 'All_stresses' ;

PlotModesStress_SUBDOMAIN(DATALOC,BasisSTRESS,MESH); 

  
end




% -------------------------------------------------
% POST-PROCESSING SELF-EQUILIBRATED MODES
% -------------------------------------------------
