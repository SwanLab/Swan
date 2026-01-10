function CONNECTbDOM = CONNECTb_dom(CONNECTb,NODESrve)

if nargin == 0
    CONNECTb = [1 2; 
                2 3 
                3  4
                4 5 
                5 2 
                5 6 
                6 1] ; 
            
            NODESrve = [2 3 5 4]' ; 
    
end
% This should be programmed in a vectorized fashion 
setBelem = ElemBndSort(CONNECTb,CNb)



% ElementsBNDdom = ones(size(CONNECTb,1)) ; %  
% for inode = 1:size(CONNECTb,2)
%     SETNODES = CONNECTb(:,inode) ;  % 
%     [ELEMLOC iA]= setdiff(SETNODES,NODESrve) ; 
%     ElementsBNDdom(iA,inode) = 0 ; 
% end
% DECIDEIN = prod(ElementsBNDdom,2) ; 
% IndElementsBND = find(DECIDEIN) ; 
% CONNECTbDOM = CONNECTb(IndElementsBND,:) ; 