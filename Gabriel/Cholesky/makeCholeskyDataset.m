function out = makeCholeskyDataset(inFile,outFile,s)
d    = load(inFile);
yF   = s.data.yFeatures;
A    = mandelMap(d.columns(yF));
n    = sqrt(size(A,1));
[diagonal,lower,names] = choleskyLayout(n);
rho  = d.data(:,d.columns == "rho");
keep = rho > 0;
M    = reshape(A*d.data(keep,yF)',n,n,[]);
if s.output.densityFactor
    M = M./reshape(rho(keep),1,1,[]);
end
Z    = zeros(size(M,3),numel(diagonal)+numel(lower));
for k = 1:size(M,3)
    [L,flag] = chol(M(:,:,k),'lower');
    if flag ~= 0
        error('makeCholeskyDataset:NotSPD','Sample %d of %s is not positive definite.',k,inFile);
    end
    Z(k,:) = [log(L(diagonal))',L(lower)'];
end
out            = d;
out.data       = d.data(keep,:);
out.data(:,yF) = Z;
out.columns    = d.columns;
out.columns(yF) = names;
out.tensors    = d.tensors(:,:,:,:,keep);
out.rhoNominal = d.rhoNominal(keep);
out.source     = d.source(keep);
save(outFile,'-struct','out');
fprintf('%s: %d of %d samples (rho > 0) -> %s\n',inFile,nnz(keep),numel(keep),outFile);
end