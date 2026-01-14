function ListGauss =  small2large(IndElements,ngaus)
nelem = length(IndElements) ; 
ListGauss = zeros(ngaus*nelem,1) ; 
for igaus = 1:ngaus
    ListGauss(igaus:ngaus:end) = ngaus*(IndElements-1) + igaus ;
end