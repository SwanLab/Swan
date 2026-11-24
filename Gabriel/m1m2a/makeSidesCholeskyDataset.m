function out = makeSidesCholeskyDataset(inFile,outFile,s)
d    = load(inFile);
xF   = s.data.xFeatures;
yF   = s.data.yFeatures;
A    = mandelMap(d.columns(yF));
n    = sqrt(size(A,1));
[diagonal,lower,names] = choleskyLayout(n);
factor = s.output.densityFactor.fun(d.data(:,xF));
keep   = factor > 0;
M      = reshape(A*d.data(keep,yF)',n,n,[]);
M      = M./reshape(factor(keep),1,1,[]);
Z      = zeros(size(M,3),numel(diagonal)+numel(lower));
for k = 1:size(M,3)
    [L,flag] = chol(M(:,:,k),'lower');
    if flag ~= 0
        error('makeSidesCholeskyDataset:NotSPD','Sample %d of %s is not positive definite.',k,inFile);
    end
    Z(k,:) = [log(L(diagonal))',L(lower)'];
end
out             = d;
out.data        = d.data(keep,:);
out.data(:,yF)  = Z;
out.columns     = d.columns;
out.columns(yF) = names;
out.rhoNominal  = d.rhoNominal(keep);
out.xi          = d.xi(keep);
out.tensors     = d.tensors(:,:,:,:,keep);
out.source      = d.source(keep);
save(outFile,'-struct','out');
fprintf('%s: %d of %d samples (density factor > 0) -> %s\n',inFile,nnz(keep),numel(keep),outFile);
end