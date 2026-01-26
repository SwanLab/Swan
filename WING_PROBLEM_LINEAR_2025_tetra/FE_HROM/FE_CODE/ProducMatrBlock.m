function CB = ProducMatrBlock(C,B)
% Given the matrices C = [C1;C2..Cn]  and B = [B1;B2 ...Bn], where
% size(Ci,1)=size(Ci,2) =m = size(Bi,1),  ProducMatrBlock returns the matrix defined
% by
% CB = [C1*B1;C2*B2 ;...; Cn*Bn]
% J.A. Hernández, jhortega@cimne.upc.edu , 27 Oct 2015
%-------------------------------------------------------

if nargin == 0
    C1 = [1 2; 3 4] ; C2 = [5 6; 7 8] ;
    C  =[C1;C2] ;
    B1 = [1 1 2 2  3  3; 0 1 0 1 0 1];
    B2 = B1.^2 ;
    B = [B1;B2];
    m = size(B1,1) ;
    CB_check= [C1*B1;C2*B2];
    load('tmp1.mat')
end

m = size(C,2) ;
n = size(C,1)/m ;
p = size(B,2) ;
nzmaxCB = nzmax(B) ;  % Maximum number of zeros
CB = zeros(n*m,p) ;
for i=1:m
    iglo = i:m:m*n ;
    for j=1:m
        jglo = j:m:m*n;
        % for k=1:p
        %    CB(iglo,k) = CB(iglo,k) + C(iglo,j).*B(jglo,k) ;
     %   CBij = 
        CB(iglo,:) = CB(iglo,:) + bsxfun(@times,B(jglo,:),C(iglo,j)) ; 
        %         nzmaxLOC = nzmax(CBij) ;
        %         CB = CB + sparse(ii,jj,ss,n*m,p,nzmaxLOC) ;
        
        % end
    end
end
%CB-CB_check
