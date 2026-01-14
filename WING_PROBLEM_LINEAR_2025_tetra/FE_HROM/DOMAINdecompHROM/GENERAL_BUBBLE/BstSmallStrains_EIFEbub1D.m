function Bst_F= BstSmallStrains_EIFEbub1D(MESH,nstrain,Bmat_all)
%  Adaptation of BstSmallStrains_EIFEbub.m  to 1D problems. In turn, this
%  is an adaption of  Bst_F= BstSmallStrains_EIFE.m to bubble modes
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/04_GeneralTheory.mlx
% JAHO, 3-Oct-2023, UPC, Terrassa/13-Oct-2023, Buenas Migas (Diag.),
% Barcelona / 13-May-2024, Honest Greens, Tuset, Barcelona 
% -----------------------------------------------------------------------------------------------------------------
if nargin ==0 
    load('tmp2.mat')
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
ngaus = size(Bmat_all{1},1)/nstrain ;
disp('Assembly of Bst (stacked B-matrix)...')
[nnode,ndim ]= size(MESH.COOR) ;
% if ndim == 2
%     nstrainF = 4;
%     IDENTITY = [1;1;0;0] ;
% elseif ndim ==3 
%     nstrainF = 9 ; 
%     IDENTITY = [1;1;1;0;0;0;0;0;0] ; 
% end

NDOFS_pernode = MESH.NDOFS_pernode ; 

Bst_F = AssemblyBGlobal(Belem_F,nstrain,nelem,nnodeE,NDOFS_pernode,ngaus,MESH.CN,nnode) ;

% ELIMINATING GHOST DOFS

Bst_F = Bst_F(:,MESH.DOFS_TO_KEEP) ; 

IDENTITY_F  =  [] ; 


