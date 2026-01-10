function [rDEF ]= AmplitudeReactions_onlybeam(DATAROM,MESH1D,a,fextBEAMr,ndim) 

if nargin ==0 
    load('tmp1.mat')
end

%V = DATAROM.BasisInt ;  % Interface modes matrix
%ndim = size(V,2) ; % Number of entries for each node
nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
nnodeE = 2; % Number of nodes per element (number of interfaces per element)

rDEF = cell(nelem,1) ;  % Amplitude self-equilibrated reaction modes 
rRB = cell(nelem,1) ;    % Amplitude resultant reaction modes 

T_1 = DATAROM.BasisInt'*DATAROM.BasisRdef(DATAROM.f1,:) ; 
T_2 = DATAROM.BasisInt'*DATAROM.BasisRdef(DATAROM.f2,:) ; 

KT{1} = DATAROM.Kbeam*T_1' ; 
KT{2} = DATAROM.Kbeam*T_2' ; 


for ielem = 1:nelem
    
    CNlocNOD = MESH1D.CN(ielem,:) ;
    % CNloc = Nod2DOF(CNlocNOD,ndim) ;
    rDEF{ielem} = -fextBEAMr{ielem}  ;
    for inode = 1:nnodeE
        NODE = CNlocNOD(inode) ;
        DOFS = Nod2DOF(NODE,ndim) ;
        rDEF{ielem} = rDEF{ielem} + KT{inode}*a(DOFS)  ; 
    end
    
end
