function T_i = RotateMatricesProduct_interface(R,BasisRdef_f,V_i)

if nargin == 0
    load('tmp4.mat')
end
if isempty(R)
    T_i = V_i'*BasisRdef_f ;
else    
    ndim = size(R,1) ;    
    BasisRdef_rot = zeros(size(BasisRdef_f)) ;
    for imodes = 1:size(BasisRdef_f,2)        
        LOCV = (R*reshape(BasisRdef_f(:,imodes),ndim,[])) ;
        BasisRdef_rot(:,imodes) = LOCV(:) ;        
    end    
    T_i = V_i'*BasisRdef_rot ;  
end