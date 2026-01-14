function [DOFr,dR,G,DOFm] = BCsCoarseScaleGeneral(ISSCCROT,ISSCC,NODES_LINES,ilines,DISPLOC,ndim,V,DATAIN,MESH2D,DATA_REFMESH)
if nargin == 0
    load('tmp.mat')
end

% Method after 16-July-2019  (only when some variable is unrestrained)
% ATTENTION ! IN THIS METHOD IS ASSUMED THAT THE ROTATIONS ARE
% UNCONSTRAINED. THAT IS, WE ONLY USE THE TRANSLATIONAL ENTRIES OF DISPLOC
if any(ISSCCROT==0)
    error('IN THIS METHOD (BoundaryConditionsUnsconstrainedDOFS) IS ASSUMED THAT THE ROTATIONS ARE  UNCONSTRAINED')
end

DOFStotal = small2large(NODES_LINES{ilines},ndim) ;  % Degrees of freedom associated to the entity under study
% Such DOFs are
% calculated using ndim
% = max (ndim 2 FACES)

nnodes = length(NODES_LINES{ilines}) ;  % Number of DOFS
% Known DOFs

ndimLOC = length(DISPLOC) ;  % Number of rigid body modes


INDX = find(ISSCC == 0) ; % Constrained RIGID BODY DOFs
% Basis matrix associated to this index
% First we have to determine which is the face associated to the line/face under
% study
[~,iface]   = find(MESH2D.CN==NODES_LINES{ilines}(1));
Vloc = V{iface} ;   % Basis matrix (interface modes )
% The prescribed displacements is therefore
DISPLOCrb = cell2mat(DISPLOC(INDX)) ;
u = Vloc(:,INDX)*DISPLOCrb;
% On the other hand, we have that
INDkeep = [] ;  % Indices associated to the restricted DOFs 
for iii = 1:length(INDX)
    INDkeep = [INDkeep,INDX(iii):DATAIN.ndimSP:size(Vloc,1)] ;
end
INDkeep = sort(INDkeep(:)) ;
Vbar = Vloc(INDkeep,:) ;
u = u(INDkeep,:) ;
% Rank of the matrix  (relatively large tolerance = 1e-3)
tol = 1e-3;
DATAIN = DefaultField(DATAIN,'BoundaryConditionsUnsconstrained_TOLERANCE',tol) ; 
tol = DATAIN.BoundaryConditionsUnsconstrained_TOLERANCE ; 


[~,DOFslave]=licols(Vbar,tol)  ;  % Slave DOFs 
DOFmast = 1:size(Vbar,2) ;
DOFmast(DOFslave) = [] ;  % Master DOFs 

Vs = Vbar(:,DOFslave) ;  % Decomposition of Vbar 
Vm = Vbar(:,DOFmast) ;


M1d = DATA_REFMESH.GeometricMassMatrixInterface{iface} ;
ndimLOC = length(INDX) ;
Mbar  = sparse(size(M1d,1)*ndimLOC,size(M1d,2)*ndimLOC) ;
for idim = 1:ndimLOC
    Mbar(idim:ndimLOC:end,idim:ndimLOC:end) = M1d ;
end

% Mbar = speye(size(Mbar)) ; 

uBARloc = (Vs'*Mbar*Vs)\(Vs'*Mbar*u) ;   % For one single coarse scale node
Gloc = (Vs'*Mbar*Vs)\(Vs'*Mbar*Vm) ;     % For one single coarse scale node 

% % Therefore  ... List of slave DOFs 
% % ---------------------------------
% if ilines==1
% warning('Temporary')
% % c = Vbar(:,6)\Vbar(:,4) ; 
% % Gloc = sparse(size(Gloc)); 
% % Gloc(4,2) = c; 
%  Gloc(1:3,2) = 0 ; 
% end

 DOFr = [] ; % DOFs with prescribed motion (subset of DOFStotal)
 DOFm = [] ; 
 dR = [] ; 
 
 nSLV = nnodes*length(DOFslave) ; 
 nMAST =  nnodes*length(DOFmast) ; 
 G = sparse(nSLV,nMAST) ; 
 
 for inode = 1:nnodes
     INDREF = (inode-1)*ndim ; 
     indDOFrLOC = DOFslave+INDREF ;      
     indDOFmLOC = DOFmast+INDREF ;       
     DOFr = [DOFr ; DOFStotal(indDOFrLOC)] ; 
     DOFm = [DOFm ; DOFStotal(indDOFmLOC)] ; 
     dR = [dR ; uBARloc] ; 
     
     IINIslave  = (inode-1)*length(DOFslave) +1 ; 
     IFINslave = inode*length(DOFslave) ; 
     
     IINImast  = (inode-1)*length(DOFmast) +1 ; 
     IFINmast = inode*length(DOFmast) ; 
     
     G(IINIslave:IFINslave,IINImast:IFINmast) = Gloc ; 
     
 end
 
 


