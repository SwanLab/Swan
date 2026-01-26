function P = AssemblyPinterfacesNON(DATAROM,MESH1D,DATAIN,FORCES,DATA_REFMESH,ndim)

% Copy of AssemblyPinterfaces. Extension to time-dependent problems 
if nargin == 0
    load('tmp2.mat')
end



nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
nnodeE = 2; % Number of nodes per element (number of interfaces per element)

% Node left/right end 
% -------------

NODE{1} = MESH1D.LEFT_END_NODE ; 
NODE{2} = MESH1D.RIGHT_END_NODE ; 
FORCESloc{1} = FORCES.LEFT_END;  
FORCESloc{2} = FORCES.RIGHT_END;  
nnodesBND = 2; 
P.VALUE = cell(1,nnodesBND) ;  
P.VALUE(:) = {zeros(nnode*ndim,1)} ; 
P.TIME_FACTOR{1} = FORCES.TIME_FACTORS.LEFT_END;  
P.TIME_FACTOR{2} = FORCES.TIME_FACTORS.RIGHT_END;  


for inode =1:nnodesBND
DOFS = small2large(NODE{inode},ndim) ; 
%nforces =length(FORCESloc{inode}) ; 
P.VALUE{inode}(DOFS ) = FORCESloc{inode} ; 

 

end