function Unew = QUprod(Q,U,alphaA)

if nargin ==0
    load('tmp1.mat')
end

if nargin == 2
   alphaA = [] ;
    
    
end



 
[alpha beta]= cellfun(@size,Q);

 if any(alpha==0)
     alpha = alphaA ; 
  end

M = sum(alpha);
Unew = zeros(M,size(U,2)) ;
iacum = 1;
iacumNEW = 1;
for i =1:length(Q)
    
    if  ~isempty(Q{i})
        nrows = size(Q{i},1) ;
        indLLLnew = iacumNEW:(nrows+iacumNEW-1);
        ncols = size(Q{i},2) ;
        indLLL = iacum:(ncols+iacum-1);
        
        Unew(indLLLnew,:) = Q{i}*U(indLLL,:) ;
        iacum = iacum+ncols ;
        
    else
        nrows =  alpha(i) ; %(Q{i},1) ;
        indLLLnew = iacumNEW:(nrows+iacumNEW-1);
    %    ncols = size(Q{i},2) ;
     %   indLLL = iacum:(ncols+iacum-1);
        
      %  Unew(indLLLnew,:) = Q{i}*U(indLLL,:) ;
      %  iacum = iacum+ncols ;
        
    end
    
    iacumNEW = iacumNEW+nrows ;
end