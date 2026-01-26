function uBAR = DomainRedDirichletConditionss(uBAR_ENDS,ALPHA_ENDS,f1,f2,BasisUrb,ndim)

  
% Dirichlet conditions left and right end
% ------------------------------
uBAR =cell(1,length(uBAR_ENDS));
for i = 1:length(uBAR_ENDS)
     if i==1
            FACE = f1 ; 
        else
            FACE = f2 ; 
        end
    if ALPHA_ENDS(i) ==0
        uBAR{i} = zeros(size(FACE))  ; 
    else
       
        uBAR_0 = uBAR_ENDS{i} ; 
        if ndim ==2
            uBAR{i} = zeros(size(FACE)) ;
          %  ndim = size(COOR,2) ;
            uBAR{i}(1:2:end) = uBAR_0(1) ;
            uBAR{i}(2:2:end) = uBAR_0(2) ;
            uBAR{i}= uBAR{i} + uBAR_0(3)*BasisUrb(FACE,3) ;
        else
            uBAR{i} = zeros(size(FACE)) ;
            
            uBAR{i}(1:3:end) = uBAR_0(1) ;
            uBAR{i}(2:3:end) = uBAR_0(2) ;
            uBAR{i}(3:3:end) = uBAR_0(3) ;
            uBAR{i}= uBAR{i} + uBAR_0(4)*BasisUrb(FACE,4) ;
            uBAR{i}= uBAR{i} + uBAR_0(5)*BasisUrb(FACE,5) ;
            uBAR{i}= uBAR{i} + uBAR_0(6)*BasisUrb(FACE,6) ;
        end
    
        
    end
end