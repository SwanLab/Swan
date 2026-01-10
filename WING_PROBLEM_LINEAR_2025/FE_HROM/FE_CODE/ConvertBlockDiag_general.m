function [Cdiag,irowsORIG_all, icolsDIAG_all]= ConvertBlockDiag_general(C,nrows,irowsORIG_all,icolsDIAG_all)

if nargin == 2 
    irowsORIG_all = [] ; icolsDIAG_all = []   ; 
end
% Given  matrices C = [C1;C2..Cn]   where
% size(Ci,1)=size(Cj,1) =m  ConvertBlockDiag_general returns (in sparse format) the matrix defined
% as
%  Cdiag = diag(C1,C2,...)
% J.A. Hernández,   30-Jun-2021 . See FLOAT_IMPLE.pdf
%-------------------------------------------------------

if nargin == 0
    C1 = [1 2 3; 4 5 6] ; C2 = [7 8 9; 10 11 12] ; C3 = [13 14 15; 16 17 18];
    C  =[C1;C2;C3] ;
    nrows = size(C1,1);
end



VECT  =1 ;

ncolsMAT = size(C,2) ;
nrowsMAT = nrows;
nmat = size(C,1)/nrowsMAT ;
n = ncolsMAT*nmat ;
m = nrowsMAT*nmat;


if   VECT ==0
    irowsORIG_all = zeros(prod(size(C)),1) ;
    icolsORIG_all = zeros(prod(size(C)),1) ;
    irowsDIAG_all = zeros(prod(size(C)),1) ;
    icolsDIAG_all = zeros(prod(size(C)),1) ;
    
    Cdiag = sparse(m,n) ;
    
    
    
    iacum = 0 ;
    for e = 1:nmat
        for irows = 1:nrowsMAT
            irowsORIG = (e-1)*nrowsMAT + irows ;
            irowsDIAG = irowsORIG ;
            for icols = 1:ncolsMAT
                icolsORIG =  icols ;
                icolsDIAG = (e-1)*ncolsMAT + icols ;
                
                %    Cdiag(irowsDIAG,icolsDIAG) = C(irowsORIG,icolsORIG) ;
                iacum = iacum+1 ;
                irowsORIG_all(iacum) = irowsORIG   ;
                icolsORIG_all(iacum) =  icolsORIG   ;
                
                irowsDIAG_all(iacum) = irowsDIAG   ;
                icolsDIAG_all(iacum) =  icolsDIAG   ;
            end
        end
    end
    C = C' ;
    C= C(:) ;
    Cdiag = sparse(irowsDIAG_all,icolsDIAG_all,C) ;
    
    Cdiag = full(Cdiag);
    
else
    
    
    if   isempty(icolsDIAG_all)
        irowsORIG_all = 1:m ;
        irowsORIG_all = repmat(irowsORIG_all,ncolsMAT,1) ;
        irowsORIG_all = irowsORIG_all(:);
        
        IND = (1:n)' ;
        IND = reshape(IND,ncolsMAT,[]) ;
        icolsDIAG_all = zeros(ncolsMAT,size(IND,2)*nrowsMAT) ;
        for icols = 1:nrowsMAT
            INI = icols:nrowsMAT:size(icolsDIAG_all,2) ;
            icolsDIAG_all(:,INI) = IND ;
        end
        icolsDIAG_all = icolsDIAG_all(:) ;
        
            
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nzmax = prod(size(C)) ;
    s =reshape(C',nzmax,1);
    %
    Cdiag= sparse(irowsORIG_all,icolsDIAG_all,s,m,n,nzmax)  ;
    
end

%full(Cdiag)