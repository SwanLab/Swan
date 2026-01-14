function  [DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = PrescFluctQ4(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC,Nshape,rnodLOC,COMB_FACES)
% JAHO, 5-Dec-2023
% See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/14_GIVEN_FLUCT_Q4.mlx
if nargin == 0
    load('tmp.mat')
end
[b,c,C,f,F,f_wcorners,F_wcorners]  =FaceCornersDOFS_Q4(rnodLOC,ndim,MESH,DATALOC,COMB_FACES) ;

INFO_PERIODIC_CONDITIONS.boundaryDOFS = b ;
INFO_PERIODIC_CONDITIONS.cornerDOFS_glo = c ;
INFO_PERIODIC_CONDITIONS.cornerDOFS_loc = C ;
INFO_PERIODIC_CONDITIONS.facesDOFS_glo = f ;
INFO_PERIODIC_CONDITIONS.facesDOFS_loc = F ;
INFO_PERIODIC_CONDITIONS.facesDOFS_glo_withCORNERS = f_wcorners ;
INFO_PERIODIC_CONDITIONS.facesDOFS_loc_withCORNERS = F_wcorners ;

% Matrix of modes
nmodes = ndim*size(Nshape,2) ;
ndofsLOC = ndim*size(Nshape,1) ;
V =zeros(ndofsLOC,nmodes) ;

for innode = 1:size(Nshape,2)
    for idim = 1:ndim
        imode = ndim*(innode-1)+idim ;
        V(idim:ndim:end,imode) = Nshape(:,innode) ;
    end
end


% \item Let $\alphaB_i \subset \{1,2 \ldots P_i \} $  be a set of $s_i$ entries of $\F_i$ such that $\U^i_{\alpha_i}$ is invertible. 
% With this decomposition at our disposal, we can split \refeq{eq:9,...lpmed} into two equations 
 U = DATALOC.Ufluct; % FLuctuations faces 1 and 2
 
alphaB = cell(size(U)) ; 
alphaB = cell(size(U)) ; 

for ifacePER =1:length(U)
   [~,alphaB{ifacePER}]=licols(U{ifacePER}') ; %
 betaB{ifacePER} = setdiff(1:size(U{ifacePER},1),alphaB{ifacePER}) ; 

end

  

% %  MASTER
% \begin{equation}
%  \d_M =   \colcuatro{\d_{f_1(\alpha_1)}}{\d_{f_2(\alpha_2)}}{\d_{f_3(\alpha_1)}}{\d_{f_4(\alpha_2)}}
% \end{equation}

M = [f{1}(alphaB{1});f{2}(alphaB{2});f{3}(alphaB{1});f{4}(alphaB{2})] ;  % Master DOFs
% %
% \item Slave DOFs
%
% \begin{equation}
%  \d_S =   \colcinco{\d_{f_1(\beta_1)}}{\d_{f_2(\beta_2)}}{\d_{f_3(\beta_1)}}{\d_{f_4(\beta_2)}}{\d_{c}}
% \end{equation}
S = [f{1}(betaB{1});f{2}(betaB{2});f{3}(betaB{1});f{4}(betaB{2});c] ;  % Master DOFs

% \begin{equation}
%  \d_M  = \G \d_S + \uBARb
% \end{equation}
%

%
% where 
%  \begin{equation}
%   \J^i \defeq \U^{i}_{\beta_i} \U^{i^{-1}}_{\alpha_i}
%  \end{equation}
J =cell(size(U)); 
for i =1:length(U)
   J{i} = U{i}(betaB{i},:)*inv(U{i}(alphaB{i},:)) ; 
end
% where
% %
% \begin{equation}
%  \G \defeq \begin{bmatrix}
%   \J^1  & \zero &  \zero &  \zero     \\
%    \zero & \J^2  & \zero &  \zero       \\
%     \zero  &  \zero & \J^1  & \zero     \\
%  \zero &  \zero  &  \zero & \J^2        \\
%   \zero  & \zero &  \zero &  \zero      \\
%  \end{bmatrix}  
% \end{equation}


%
G = sparse(length(S),length(M)) ;
% Block G_11
iniROWS = 1; 
finROWS = length(alphaB{1}) ; 
iniCOLS = 1; 
finCOLS = length(betaB{1}) ; 
G(iniCOLS:finCOLS,iniROWS:finROWS) = J{1}  ; 
% Block G_22
iniROWS = finROWS+1; 
finROWS = iniROWS+length(alphaB{2})-1 ; 
iniCOLS = finCOLS+1; 
finCOLS = iniCOLS + length(betaB{2})-1 ; 
G(iniCOLS:finCOLS,iniROWS:finROWS) = J{2}  ; 
% Block G_33
iniROWS = finROWS+1; 
finROWS = iniROWS+length(alphaB{1})-1 ; 
iniCOLS = finCOLS+1; 
finCOLS = iniCOLS + length(betaB{1})-1 ; 
G(iniCOLS:finCOLS,iniROWS:finROWS)  = J{1}  ; 
% Block G_44
iniROWS = finROWS+1; 
finROWS = iniROWS+length(alphaB{2})-1 ; 
iniCOLS = finCOLS+1; 
finCOLS = iniCOLS + length(betaB{2})-1 ; 
G(iniCOLS:finCOLS,iniROWS:finROWS) = J{2}  ; 

% %%%%%%%%%%%%%%%%%%%%
%   \begin{equation}
%   \b^i  =  \V_{F_i(\beta_i)} \a  -  \U^{i}_{\beta_i}  \U^{i^{-1}}_{\alpha_i} \V_{F_i(\alpha_i)} \a
%  \end{equation}
a = DIRICHLET.AMPLITUDE;
 


b_1 = V(F{1}(betaB{1}),:)*a  - U{1}(betaB{1},:)*(U{1}(alphaB{1},:)\(V(F{1}(alphaB{1}),:)*a))   ; 
b_2 = V(F{2}(betaB{2}),:)*a  - U{2}(betaB{2},:)*(U{2}(alphaB{2},:)\(V(F{2}(alphaB{2}),:)*a))   ;
b_3 = V(F{3}(betaB{1}),:)*a  - U{1}(betaB{1},:)*(U{1}(alphaB{1},:)\(V(F{3}(alphaB{1}),:)*a))   ; 
b_4 = V(F{4}(betaB{2}),:)*a  - U{2}(betaB{2},:)*(U{2}(alphaB{2},:)\(V(F{4}(alphaB{2}),:)*a))   ;

% 
% \begin{equation}
%  \uBARb \defeq \colcinco{\b^1}{\b^2}{\b^3}{\b^4}{\V_C \a}
% \end{equation}

uBARb = [ b_1;
    b_2; 
    b_3;
    b_4;
    V(C,:)*a] ;

INFO_PERIODIC_CONDITIONS.ImposedDisplacement = V*a ;




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
DISP_CONDITIONS.G = G ; %
DISP_CONDITIONS.DOFr = S ; % SLAVE DOFS
DISP_CONDITIONS.DOFm = M ; % MASTER DOFS
DISP_CONDITIONS.DOFl = L ; % Unknown DOFS

dR.U = uBARb ;
dR.a = DIRICHLET.TIMEFUN(DATA.STEPS) ;
DISP_CONDITIONS.dR = dR ;