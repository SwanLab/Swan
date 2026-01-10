function gaussFIB = GaussPointsAssociatedElements(ngausE,elemFIB)
gaussFIB = ngausE*repmat(elemFIB',ngausE,1) ; 
SUMind = (ngausE-1):-1:0 ; 
gaussFIB = bsxfun(@plus,gaussFIB,-SUMind');
gaussFIB = gaussFIB(:) ; 