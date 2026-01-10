function D_QrotINI_byROWS = Qrot_by_rows_element(Qrot,INDEXsparseROTmat)
% Construction of rotation matrices per rows (vectorized implementation)
% JAHO, 25-Feb-2025, Tuesday, 14:15 UPC, Terrassa .
% % /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/06_COROT_3D.mlx
%  \DiagC{[\Qrot]_i} = \diagOL{[\QrotE{1}]_i}{[\QrotE{2}]_i}{\cdots}{[\QrotE{e}]_i} ,  \hspace{0.5cm}  i  =1,2,3

if nargin == 0
    load('tmp1.mat')
    INDEXsparseROTmat = [] ;
    Qrot = [Qrot;Qrot] ;
elseif nargin == 1
    INDEXsparseROTmat = [] ;
end
ndim = 3;
D_QrotINI_byROWS = cell(1,ndim) ;

for idim = 1:ndim
    D_QrotINI_byROWS{idim} = Qrot_by_rows_element_idim(Qrot,idim,INDEXsparseROTmat)  ;
end




end
 
function D_QrotINI_byROWS = Qrot_by_rows_element_idim(Qrot,idim,INDEXsparseROTmat)
ndim = 3;
rows_selected = idim:ndim:size(Qrot);
Qrot_row = Qrot(rows_selected,:) ;
s_Qrot = Qrot_row(:);
nelem  =size(Qrot_row,1) ;
nzmax = length(s_Qrot) ;

if isempty(INDEXsparseROTmat)
    
    [irows,icols] = IndexesQrotROWS(nelem,ndim) ;
    D_QrotINI_byROWS = sparse(irows,icols,s_Qrot,nelem,nelem*ndim,nzmax) ;
else
    D_QrotINI_byROWS = sparse(INDEXsparseROTmat.IROWS,INDEXsparseROTmat.ICOLS,s_Qrot,nelem,nelem*ndim,nzmax) ;
end



end
