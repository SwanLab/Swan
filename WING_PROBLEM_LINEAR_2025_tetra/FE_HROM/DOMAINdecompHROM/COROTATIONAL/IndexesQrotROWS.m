function  [irows,icols] = IndexesQrotROWS(nelem,ndim)
 irows  = repmat(1:nelem,1,ndim) ;
irows = irows(:) ;
icols  = repmat(1:ndim,nelem,1) ;
addROWS = [0:ndim:(nelem-1)*ndim]' ;
icols = bsxfun(@plus,icols,addROWS) ;
icols = icols(:) ;