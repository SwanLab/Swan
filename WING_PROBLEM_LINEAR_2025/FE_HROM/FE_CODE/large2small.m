function ListElements =  large2small(IndGauss,ngaus)

if nargin == 0
    ngaus = 8 ; 
    IndGauss =[1 2 3 4 9 20 ] ; 
    
   
end

DIV = IndGauss/ngaus ;
ListElements = ceil(DIV) ; 

ListElements = unique(ListElements) ; 

 if size(ListElements,1) ==1
    ListElements = ListElements'; 
end


% nelem = length(IndElements) ; 
% ListGauss = zeros(ngaus*nelem,1) ; 
% for igaus = 1:ngaus
%     ListGauss(igaus:ngaus:end) = ngaus*(IndElements-1) + igaus ;
% end