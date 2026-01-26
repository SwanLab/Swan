function CB = ProducMatrBlockGen(C,B,nblock)
% Given the matrices C = [C1;C2..Cn]  and B = [B1;B2 ...Bn], where
%  size(Ci,2) =m = size(Bi,1),  ProducMatrBlockGen returns the matrix defined
% by
% CB = [C1*B1;C2*B2 ;...; Cn*Bn]
% J.A. Hernández, jhortega@cimne.upc.edu , 28 Oct 2015
%-------------------------------------------------------

if nargin == 0
    C1 = [ 1 2; 3 4; -3 -4] ; C2 = [ 5 6; 7 8; -7 -8] ;
    C  =[C1;C2] ;
    B1 = [1 1 2 2  3  3; 0 1 0 1 0 1];
    B2 = B1.^2 ;
    B = [B1;B2];
    m = size(B1,1) ;
    CB_check= [C1*B1;C2*B2];
    nblock = 2 ;
    %   load('tmp1.mat')
end

my = size(C,2) ;
n = nblock ;
mx = size(C,1)/n ;
p = size(B,2) ;
CB = zeros(n*mx,p) ;
for i=1:mx
    iglo = i:mx:mx*n ;
    for j=1:my
        jglo = j:my:my*n;
        % for k=1:p
        %    CB(iglo,k) = CB(iglo,k) + C(iglo,j).*B(jglo,k) ;
        CBij = bsxfun(@times,B(jglo,:),C(iglo,j)) ;
        CB(iglo,:) = CB(iglo,:) +  CBij ;
        %         nzmaxLOC = nzmax(CBij) ;
        %         CB = CB + sparse(ii,jj,ss,n*m,p,nzmaxLOC) ;
        
        % end
    end
end
%CB-CB_check
