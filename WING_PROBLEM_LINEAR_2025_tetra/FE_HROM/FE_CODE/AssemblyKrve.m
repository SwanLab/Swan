function [K,DOFsKEEP] = AssemblyKrve(DATAROM,MESH2D,DATAIN,ndim)
% JAHO, 4-July-2018/24-July-2018
% Reduced-order Stiffness matrix
% All interfaces are assumed to have the same number of modes 

if nargin == 0
    load('tmp1.mat')
end

 
 



nnodeE = length(ndim) ; 
DOFsKEEP = [] ; 
if any(ndim-ndim(1))
    [K,DOFsKEEP] = AssemblyKrve_notequal(DATAROM,MESH2D,DATAIN,ndim) ; 
else
   % Standard assembly all interfaces posses the same number of modes
    [K] = AssemblyKrve_ndimeq(DATAROM,MESH2D,DATAIN,ndim(1),nnodeE) ; 
end
 