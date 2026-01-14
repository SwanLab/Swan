function  [DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] =...
    Periodic_1directpz(DIRICHLET,DATA,ndim,MESH,DATALOC)
% Goal. Determine DOFr and    dR(t) .
% Periodicity conditions in one direction
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/01_meta_1D.mlx
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
% 
% 
%
% JAHO- 14,31-May-2O24, Terrassa,  UPC
% ------------------------------------------------------
if nargin == 0
    load('tmp3.mat')
end
INFO_PERIODIC_CONDITIONS = [] ;

NODES_FPER = cell(1,2) ;
[NODES_FPER{1},NODES_FPER{2}] = SetNodesPeriodic1D(DIRICHLET,MESH) ;


f = cell(1,2) ;
c = cell(size(f)) ;

for indLOC = 1:length(DIRICHLET.FacesPeriodicity1D)
    nodes_FACE = NODES_FPER{indLOC}  ;
    % Corner 1
    nc= DIRICHLET.NodesCornersPeriodicity1D(indLOC) ;
    nodes_FACE_wc  = setdiff(nodes_FACE,nc,'stable') ;
    if length(nodes_FACE_wc) == length(nodes_FACE)
        error(['Node ',num2str(nc),' does not belong to face ',num2str()])
    end
    f{indLOC} = small2large(nodes_FACE_wc,ndim ) ;
    c{indLOC} = small2large(nc,ndim ) ;
    
end


% REMAINING DOFs to be set to zero 
% 
for iface = 1:length(DIRICHLET.FacesZeroDOFs)
    indFACEloc = DIRICHLET.FacesZeroDOFs(iface);  
       NODES_loc = MESH.NODES_FACES{indFACEloc} ;
       All_DOFS_loc = small2large(NODES_loc,ndim )  ; 
       SELECTED_DOFS_zer_constraint = [];
       indDOFsLOC = DIRICHLET.DOFs_loc_constrained{iface}; 
       for iii = 1:length(indDOFsLOC)
           iLOC = indDOFsLOC(iii) ; 
           SELECTED_DOFS_zer_constraint = [SELECTED_DOFS_zer_constraint; All_DOFS_loc(iLOC:ndim:end)] ;
       end
       
end 


  %
% \begin{equation}
%  \d_M =  \d_{f{1}}
% \end{equation}
%
% \item Slave DOFs
%
% \begin{equation}
%  \d_S = \coltres{\d_{f{2}}}{d_{c{1}}}{d_{c{2}}}
% \end{equation}


M = [f{1}] ;  % Master DOFs
S = [f{2};c{1}; c{2};SELECTED_DOFS_zer_constraint] ;  % Slave DOFs

%
% \begin{equation}
%   \G \defeq  \begin{bmatrix}
%                                           \ident  \\
%                                           \zero  \\
%                                           \zero
%                                          \end{bmatrix}
%                                          \end{equation}
%
% and
%
% \begin{equation}
%  \uBARb \defeq  \coltres{\ones}{0}{\coldos{1}{0}} a  
% \end{equation}

%
G = sparse(length(S),length(M)) ;
nrowsf12 = length(f{1}) ;% + length(f{2}) ;
n_f12 = 1:nrowsf12 ;
G(n_f12,n_f12) = speye(nrowsf12,nrowsf12) ;

a = DIRICHLET.PRESCRIBED_DISP.AMPLITUDE;

MODE = zeros(size(f{1})) ; 
MODE(1:2:end) = 1; 


uBARb = [MODE
    0
    0
    1
    0
    zeros(size(SELECTED_DOFS_zer_constraint))]*a ;


nDOFS = size(MESH.COOR,1)*size(MESH.COOR,2);
I = setdiff(1:nDOFS,[M;S])' ;   % Interior DOFs

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
DISP_CONDITIONS.G = G ; %
DISP_CONDITIONS.DOFr = S ; % SLAVE DOFS
DISP_CONDITIONS.DOFm = M ; % MASTER DOFS
DISP_CONDITIONS.DOFl = L ; % Unknown DOFS

dR.U = uBARb ;
dR.a = DIRICHLET.PRESCRIBED_DISP.TIMEFUN(DATA.STEPS) ;
DISP_CONDITIONS.dR = dR ;