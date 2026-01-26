function ListGauss =  Elements2Gauss(ngaus,nelem,IndElements)

ListGauss = zeros(ngaus*nelem,1) ; 
for igaus = 1:ngaus
    ListGauss(igaus:ngaus:end) = ngaus*(IndElements-1) + igaus ;
end