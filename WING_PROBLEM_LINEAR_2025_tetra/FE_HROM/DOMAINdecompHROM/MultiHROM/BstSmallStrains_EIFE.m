function Bst_F= BstSmallStrains_EIFE(MESH,nstrain,Bmat_all)
% Adaptation BstSmallStrains, EIFE elements
% JAHO, 20-March-2023, Teknon clinic, Barcelona
if nargin ==0 
    load('tmp.mat')
end
% 
% DATALOC.DEFORMATION_GRADIENT  = 0 ; 
% MESH = DefaultField(MESH,'posgp_given',[]) ; 
% MESH = DefaultField(MESH,'weights_given',[]) ; 
% 
% DATALOC.posgp_given = MESH.posgp_given ; 
% DATALOC.weights_given = MESH.weights_given ; 
% DATALOC.StrainStressWith4Components = MESH.DATA.StrainStressWith4Components ; 

[aa,bb] = cellfun(@size,Bmat_all,'UniformOutput',false);
ngaussall = cell2mat(aa) ; 

if length(unique(ngaussall)) ~=1
    error('Option not implemented yet (all elements should have the same number of points)')
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

Bst_F = AssemblyBGlobal(Belem_F,nstrain,nelem,nnodeE,ndim,ngaus,MESH.CN,nnode) ;


IDENTITY_F  =  [] ; 


