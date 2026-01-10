function  [DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = PeriodicQ4_NOcorners...
    (DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC,Nshape,rnodLOC,COMB_FACES)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/13_Q4_nocorners.mlx
%  JAHO, 28-Nov-2023, Terrassa, UPC. 
if nargin == 0
    load('tmp.mat')
end

[b,f,F,alphaB,betaB,Z]   =FaceNOCornersDOFS_Q4(rnodLOC,ndim,MESH,DATALOC,COMB_FACES,GEOproperties) ;  

INFO_PERIODIC_CONDITIONS.boundaryDOFS = b ; 
INFO_PERIODIC_CONDITIONS.alphaB_F3 = alphaB ; 
INFO_PERIODIC_CONDITIONS.betaB_F3 = betaB ; 
 INFO_PERIODIC_CONDITIONS.facesDOFS_glo = f ; 
INFO_PERIODIC_CONDITIONS.facesDOFS_loc = F ;  


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

% 
% 
% \item \textbf{Master and slave DOFs}.  We choose as master DOFs 
% 
% \begin{equation}
%  \d_M = \coldos{\d_{f_3(\beta)} }{\d_{f_4}}
% \end{equation}
 M = [f{3}(betaB);f{4}] ;  % Master DOFs
% and as slave DOFs 
% 
% \begin{equation}
%  \d_S =  \colcuatro{ \d_{f_1(\alpha)}}{ \d_{f_1(\beta)} }{  \d_{f_3(\alpha)}}{\d_{f_2}}
% \end{equation}
S = [f{1}(alphaB);f{1}(betaB);f{3}(alphaB); f{2} ] ;
 

% --------------------------------------------------------------------------------------------------------------
%  \begin{equation}
%   \H_3 = -\Par{\Z_{F_1(\alpha)}^T  + \Z_{F_3(\alpha)}^T  }^{-1} \Par{\Z_{F_1(\beta)}^T +  \Z_{F_3(\beta)}^T}  
%  \end{equation}
F_1a = F{1}(alphaB) ;  F_3a = F{3}(alphaB) ;
F_1b = F{1}(betaB) ;  F_3b = F{3}(betaB) ;
ZF13 = -(Z(F_1a,:)' + Z(F_3a,:)');
H_3 =  ZF13\(Z(F_1b,:)' + Z(F_3b,:)');
%  \begin{equation}
%   \H_4 = -\Par{\Z_{F_1(\alpha)}^T  + \Z_{F_3(\alpha)}^T  }^{-1}  \Par{ \Z_{F_2}^T  + \Z_{F_4}^T } 
%  \end{equation}
H_4 =ZF13\(Z(F{2},:)' + Z(F{4},:)') ;  

%%%
%  -------------------------------------------------------------------------------------------------------------
%  \begin{equation}
%   \G \defeq \begin{bmatrix}
%                                                                                       \H_3 & \H_4 \\ 
%                                                                                       \ident & \zero \\
%                                                                                        \H_3 & \H_4 \\ 
%                                                                                        \zero & \ident 
%                                                                                     \end{bmatrix}  
%  \end{equation}

%%%%% 
G = sparse(length(S),length(M)) ; 

% G_11
iniROW = 1; 
finROW = size(H_3,1) ;
iniCOL = 1; 
finCOL = size(H_3,2) ; 
G(iniROW:finROW,iniCOL:finCOL) = H_3 ; 
% G_12
iniCOL = finCOL+1; 
finCOL = size(G,2) ; 
G(iniROW:finROW,iniCOL:finCOL) = H_4 ;
% G_21
iniROW = finROW + 1; 
finROW = iniROW + length(betaB)-1 ;
iniCOL = 1; 
finCOL = size(H_3,2) ; 
G(iniROW:finROW,iniCOL:finCOL) = speye(length(betaB),length(betaB)) ;

% G_31 
iniROW = finROW+1; 
finROW = iniROW + size(H_3,1)-1 ;
iniCOL = 1; 
finCOL = size(H_3,2) ; 
G(iniROW:finROW,iniCOL:finCOL) = H_3 ;

% G_32 
iniCOL = finCOL+1; 
finCOL = size(G,2) ; 
G(iniROW:finROW,iniCOL:finCOL) = H_4 ;

% G_42 
iniROW = finROW+1; 
finROW = size(G,1); 
G(iniROW:finROW,iniCOL:finCOL) = speye(length(f{2}),length(f{2})) ; ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  \begin{equation}
%   \q \defeq \Par{ \Z^T \V  -  \Z_{F_1(\alpha)}^T\Par{\V_{F_1(\alpha)} - \V_{F_3(\alpha)}} -  \Z_{F_1(\beta)}^T   \Par{\V_{F_1(\beta)} - \V_{F_3(\beta)}}  -  \Z_{F_2}^T  \Par{\V_{F_2} - \V_{F_4}}  }\a 
%  \end{equation}
a = DIRICHLET.AMPLITUDE;  
q = Z'*V - Z(F_1a,:)'*(V(F_1a,:) - V(F_3a,:))  - Z(F_1b,:)'*(V(F_1b,:)-V(F_3b,:)) -Z(F{2},:)'*(V(F{2},:)-V(F{4},:)) ;  
q = q*a ; 
%
%  \begin{equation}
%   \p \defeq \Par{\Z_{F_1(\alpha)}^T  + \Z_{F_3(\alpha)}^T  }^{-1}  \q
%  \end{equation}
p = (Z(F_1a,:)' + Z(F_3a,:)')\q; 

% %
%  \begin{equation}
%   \uBARb \defeq  \colcuatro{ \p  +  \Par{\V_{F_1(\alpha)} - \V_{F_3(\alpha)}} \a}
%          {\Par{\V_{F_1(\beta)} - \V_{F_3(\beta)}} \a}
%          {\p}
%          {\Par{\V_{F_2} - \V_{F_4}} \a}
%  \end{equation}


uBARb = sparse(size(G,1),1) ; 

% u_1
iniROW = 1; 
finROW = size(H_3,1) ;
uBARb(iniROW:finROW) = p + (V(F_1a,:) -V(F_3a,:))*a ; 
 
% u_2
iniROW = finROW + 1; 
finROW = iniROW + length(betaB)-1 ;
 uBARb(iniROW:finROW) =  (V(F_1b,:) -V(F_3b,:))*a ; 

 
% u_3
iniROW = finROW+1; 
finROW = iniROW + size(H_3,1)-1 ;
 uBARb(iniROW:finROW) =  p ; 

 

% u_4
iniROW = finROW+1; 
finROW = size(G,1); 
 uBARb(iniROW:finROW) =  (V(F{2},:) -V(F{4},:))*a ; 



     
 nDOFS = size(MESH.COOR,1)*size(MESH.COOR,2); 
 I = setdiff(1:nDOFS,[M;S])' ;   % Interior DOFs 
 
 %  \Abub \defeq \begin{bmatrix} 
%                                                         \ident &  \zero \\
%                                                         \zero & \ident \\
%                                                         \zero &  \G 
%                                                        \end{bmatrix} 
% \end{equation}




% 
% \begin{equation}
%  \uBARb \defeq  \coltres{\Par{\V_{F_1} - \V_{F_3}}  }{\Par{\V_{F_2} - \V_{F_4}}  }{\V_C} \a
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