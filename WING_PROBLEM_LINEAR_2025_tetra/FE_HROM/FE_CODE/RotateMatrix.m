function Ur = RotateMatrix(R,U)
% Rotation of the columns of matrix U (R --> Rotation matrix)

if nargin == 0
    load('tmp4.mat')
end
if isempty(R)
    Ur = U;
else
    ndim = size(R,1) ;
    Ur = zeros(size(U)) ;
    for imodes = 1:size(U,2)        
        LOCV = (R*reshape(U(:,imodes),ndim,[])) ;
        Ur(:,imodes) = LOCV(:) ;        
    end  
    
end