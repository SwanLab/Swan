function  [DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = PeriodicQ4_corners(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC,Nshape,rnodLOC,COMB_FACES)

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
%  \d_M  = \G \d_S + \uBARb
% \end{equation}
%
% where
%
% \begin{equation}
%   \G \defeq  \begin{bmatrix}
%                                           \ident &  \zero  \\
%                                           \zero &  \ident \\
%                                           \zero & \zero
%                                          \end{bmatrix}
%                                          \end{equation}
%
% and
%
% \begin{equation}
%  \uBARb \defeq  \coltres{\Par{\V_{F_1} - \V_{F_3}}  }{\Par{\V_{F_2} - \V_{F_4}}  }{\V_C} \a
% \end{equation}
%
G = sparse(length(S),length(M)) ;
nrowsf12 = length(f{1}) + length(f{2}) ;
n_f12 = 1:nrowsf12 ;
G(n_f12,n_f12) = speye(nrowsf12,nrowsf12) ;

a = DIRICHLET.AMPLITUDE;

INFO_PERIODIC_CONDITIONS.ImposedDisplacement = V*a ;

uBARb = [V(F{1},:)-V(F{3},:)
    V(F{2},:)-V(F{4},:)
    V(C,:)]*a ;


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