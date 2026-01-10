function Cglo = DefineElastMatGLO_nw(celasglo,ngaus)
% G iven the array celasglo of elasticity matrices
% celasglo(:,:,1), celasglo(:,:,2) ... celasglo(:,:,nelem),
% this functions creates the matrix 
% Cglo = [celasglo(:,:,1)  % |
%          celasglo(:,:,1) % |> ngaus times
%         ....             % |
%         celasglo(:,:,1)  % |
%         celasglo(:,:,2)  
%         celasglo(:,:,2)
%         ......
%           celasglo(:,:,2) ... 
%         ....    ]
% J.A. Hernández, jhortega@cimne.upc.edu , 27 Oct 2015
%  Modification: 8-Dec-2020: 
%-------------------------------------------------------

if nargin == 0
    celasglo = zeros(2,2,2) ; 
    A = [1 2; 3 4] ; B = [5 6; 7 8] ; 
    celasglo(:,:,1)=A ;  celasglo(:,:,2)=B ;
    ngaus = 3 ;
    wST = 2*ones(3*2,1) ; 
end


 
nelem = size(celasglo,3) ; 
nstrain = size(celasglo,1) ; 

 Cglo=zeros(nstrain*nelem*ngaus,nstrain) ; 
for istrain=1:nstrain
    for jstrain=1:nstrain        
        colIJ =  (squeeze(celasglo(istrain,jstrain,:)))';
        colIJrep = repmat(colIJ,ngaus,1); 
        colIJ = reshape(colIJrep,[],1);
        indIJ = istrain:nstrain:nstrain*ngaus*nelem ; 
        Cglo(indIJ,jstrain) =  colIJ ; 
    end
end


