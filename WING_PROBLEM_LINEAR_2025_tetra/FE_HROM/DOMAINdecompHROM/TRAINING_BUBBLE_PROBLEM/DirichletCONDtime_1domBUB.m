function [DISP_CONDITIONS] =  DirichletCONDtime_1domBUB(DIRICHLET,DATA,ndim,MESH,MODES) ; 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/08_IndividualTRAINbub.mlx
% Boundary conditions individual domains, bubble conditions
% /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_bub.pdf
% JAHO, 29-Oct-2023, Barcelona, Balmes 185. 
% MODIFED ON NOV-19 2023. Periodic boundary conditiosn 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/10_PeriodicQUADL.mlx

if nargin == 0
    load('tmp.mat')
end


PERIODIC = 1; 

if PERIODIC == 1 
% PERIODIC BOUNDARY CONDITIONS FOR THE BUBBLE DISPLACEMENT
% ------------------------------------------------------------------------------
% New version, 19-Nov-2023 
% % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/10_PeriodicQUADL.mlx
% Patterned after /InputDataFunctions/PeriodicQ4.m
% ---------------------------------------------------------------------------------
%b = MESH.INFO_PERIODIC_CONDITIONS.boundaryDOFS  ; 
c = MESH.INFO_PERIODIC_CONDITIONS.cornerDOFS_glo   ; 
%C = MESH.INFO_PERIODIC_CONDITIONS.cornerDOFS_loc   ; 
f = MESH.INFO_PERIODIC_CONDITIONS.facesDOFS_glo   ; 
%F = MESH.INFO_PERIODIC_CONDITIONS.facesDOFS_loc   ; 

%  MASTER
% \begin{equation}
%  \d_M = \coldos{\d_{f_3} }{\d_{f_4}}
% \end{equation}
% 
% \item Slave DOFs 
% 
% \begin{equation}
%  \d_S = \coltres{\d_{f_1}}{\d_{f_2}}{\d_c}
% \end{equation}

M = [f{3};f{4}] ;  % Master DOFs
S = [f{1};f{2}; c] ;  % Slave DOFs 

% \begin{equation}
%   \G \defeq  \begin{bmatrix}
%                                           \ident &  \zero  \\
%                                           \zero &  \ident \\
%                                           \zero & \zero
%                                          \end{bmatrix}
%                                          \end{equation}
% 
% and 
G = sparse(length(S),length(M)) ; 
nrowsf12 = length(f{1}) + length(f{2}) ; 
n_f12 = 1:nrowsf12 ; 
G(n_f12,n_f12) = speye(nrowsf12,nrowsf12) ; 

 nDOFS = size(MESH.COOR,1)*size(MESH.COOR,2); 
 I = setdiff(1:nDOFS,[M;S])' ;   % Interior DOFs 
 
 %  \Abub \defeq \begin{bmatrix} 
%                                                         \ident &  \zero \\
%                                                         \zero & \ident \\
%                                                         \zero &  \G 
%                                                        \end{bmatrix} 
% \end{equation}

L = [I; M] ; % Unknown DOFs 

A = sparse(nDOFS,length(L)) ; 

 % A_11 
iniCOL = 1; 
finCOL = length(I) ; 
A(I,iniCOL:finCOL) = speye(length(I),length(I)) ;
 % A_22 
iniCOL = finCOL+1; 
finCOL = size(A,2) ; 
A(M,iniCOL:finCOL) = speye(length(M),length(M)) ; 
% A_23
A(S,iniCOL:finCOL) = G;  

DISP_CONDITIONS.A = A ; %
DISP_CONDITIONS.J = G ; %
DISP_CONDITIONS.DOFr = S ; % SLAVE DOFS
DISP_CONDITIONS.DOFm = M ; % MASTER DOFS
DISP_CONDITIONS.DOFl = L ; % Unknown DOFS

DISP_CONDITIONS.PhiDEF = MODES.PhiDEF ; 
DISP_CONDITIONS.COARSE_SCALE_HISTORY = DIRICHLET.COARSE_SCALE_HISTORY ; 



else

%  
% PsiALL 
% \begin{equation}
%  \PsiALL \defeq  \rowdos{\PsiRES}{\PsiSE}. 
% \end{equation}
PsiALL = [MODES.PsiRBf,MODES.PsiDEFf] ; % Matrix of interface forces (resultants + self-equilibrared)
%  \textbf{Master and slave DOFs}.   Let us decompose the set of interface DOFs $\bDOFS$  
% into two disjoint subsets $\hDOFS$ and $\gDOFS$ so that $\PsiALL(\hDOFS,:)$ is square and invertible.

% Patterned after :  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/BCs_BEAMS_PERIODIC.m
% 
% DATA.ChoiceDOFs_SLAVE_RANDOM =1 ; 
% 
% if DATA.ChoiceDOFs_SLAVE_RANDOM == 1 
%     minS = 0 ; 
%     while minS < 1e-4
%       hDOFS =   randperm(size(PsiALL,1),size(PsiALL,2)) ;  
%       SS = svd(PsiALL(hDOFS,:)) ; 
%       minS = SS(end)/SS(1) ; 
%     end
% else
[~,hDOFS]=licols(PsiALL') ; %
%end
% These are the slaves DOFs, in "local" notation. In global notation
bDOFS = MESH.faceDOFSall ;
%DOFr = bDOFS(hDOFS) ; % Global notation
ListNodesSlave = large2small(bDOFS(hDOFS),ndim) ; 
disp(['List slave NODES =' ]) ; 
disp( num2str(ListNodesSlave')) ; 
gDOFS = 1:length(bDOFS) ; 
gDOFS = setdiff(gDOFS,hDOFS) ; 
% MATRIX Jhg
% In doing so, we can write the boundary condition for the bubble displacements as 
% \begin{equation}
%  \PsiALL(\hDOFS,:)^T  \dFbubBh +    \PsiALL(\gDOFS,:)^T \dFbubBg  = \zero 
% \end{equation}
% Thus, 
% \begin{equation}
% \label{eq:atdisp}
%      \dFbubBh = \Jhg \dFbubBg
% \end{equation}
% where 
% \begin{equation}
%  \Jhg \defeq - \PsiALL(\hDOFS,:)^{-T}  \PsiALL(\gDOFS,:)^T
% \end{equation}
Jhg =-PsiALL(hDOFS,:)'\PsiALL(gDOFS,:)'; 

% -------------------------
% Matrix  Abub
% -------------------------
% \begin{equation}
% \label{eq:ancel}
%  \dFbub  =  \coltres{\dFbubS }{\dFbubBg}{\dFbubBh} = \begin{bmatrix} 
%                                                         \ident &  \zero \\
%                                                         \zero & \ident \\
%                                                         \zero &  \Jhg 
%                                                        \end{bmatrix}  \coldos{\dFbubS}{\dFbubBg}
% \end{equation}
% or more compactly 
% \begin{equation}
%  \dFbub = \Abub \dFbubL
% \end{equation}
% where 
% \begin{equation}
%  \Abub \defeq \begin{bmatrix} 
%                                                         \ident &  \zero \\
%                                                         \zero & \ident \\
%                                                         \zero &  \Jhg 
%                                                        \end{bmatrix} 
% \end{equation}



nDOFS = size(MESH.COOR,1)*size(MESH.COOR,2) ; 
% Non-interface DOFs 
sDOFS = setdiff((1:nDOFS)',bDOFS) ; 
% Unknown DOFs 
lDOFS = [sDOFS;bDOFS(gDOFS)] ; 

% Initialization (sparse)
Abub = sparse(nDOFS,length(lDOFS)) ;

% Block  sDOFS, 1
iini = 1; 
ifin = length(sDOFS) ; 
Abub(sDOFS,iini:ifin) = speye(length(sDOFS))  ; 

% Block  bDOFS(gDOFS), 2
iini = length(sDOFS) +1 ; 
ifin = size(Abub,2) ; 
Abub(bDOFS(gDOFS),iini:ifin) = speye(length(gDOFS))  ; 

% Block  bDOFS(hDOFS), 2
iini = length(sDOFS) +1 ; 
ifin = size(Abub,2) ; 
Abub(bDOFS(hDOFS),iini:ifin) = Jhg  ; 

% Coarse-scale strain history 
DISP_CONDITIONS.J = Jhg ; %  
DISP_CONDITIONS.A = Abub ; %  
DISP_CONDITIONS.DOFr = bDOFS(hDOFS) ; % SLAVE DOFS 
DISP_CONDITIONS.DOFm = bDOFS(gDOFS) ; % MASTER DOFS 
DISP_CONDITIONS.DOFl = lDOFS ; % remaining DOFS 
DISP_CONDITIONS.PhiDEF = MODES.PhiDEF ; 
DISP_CONDITIONS.COARSE_SCALE_HISTORY = DIRICHLET.COARSE_SCALE_HISTORY ; 

end
 