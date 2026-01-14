function [Bst_F,IDENTITY_F]= BstLargeStrains_EIFEbub(MESH,nstrain,Bmat_all,DATA)
% Adaption  Bst_F= BstLargeStrains_EIFEbub.m to bubble modes
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/03_EIFEM_large.mlx
% JAHO, 21-May-2204, Tuesday, UPC, Terrassa.  
%---------------------------------------------------------------------
if nargin ==0 
    load('tmp3.mat')
end


[aa,bb] = cellfun(@size,Bmat_all,'UniformOutput',false);
ngaussall = cell2mat(aa) ; 

if length(unique(ngaussall)) ~=1
    error('Option not implemented yet (all elements should have the same number of Gauss points/DOFs)')
end



%[ Belem_F, wSTs, wST, XeALL, posgp] = ComputeBelemALL(MESH.COOR,MESH.CN,MESH.TypeElement,nstrain,DATALOC) ;
Belem_F = cell2mat(Bmat_all) ; 
% Computing the deformation gradient as a vector 
% F = J + I , where J = Bst_F*d
[nelem,nnodeE ]= size(MESH.CN) ; 
disp('Assembly of Bst (stacked B-matrix)...')
[nnode,ndim ]= size(MESH.COOR) ;
nstrain_F = DATA.MESH.nstrain_F ; 
ngaus = size(Bmat_all{1},1)/nstrain_F ;

 if nstrain_F == 4
   %  nstrainF = 4;
     IDENTITY = [1;1;0;0] ;
 elseif nstrain_F ==9
    % nstrainF = 9 ; 
     IDENTITY = [1;1;1;0;0;0;0;0;0] ; 
 end

NDOFS_pernode = MESH.NDOFS_pernode ; 

Bst_F = AssemblyBGlobal(Belem_F,nstrain_F,nelem,nnodeE,NDOFS_pernode,ngaus,MESH.CN,nnode) ;

% ELIMINATING GHOST DOFS

Bst_F = Bst_F(:,MESH.DOFS_TO_KEEP) ; 

 IDENTITY_F  = repmat(IDENTITY,1,nelem*ngaus) ; 
IDENTITY_F = IDENTITY_F(:) ; 


