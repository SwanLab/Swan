function out = nOrderTranspose(in)
    nDims = 1:length(size(in));
    nDims(1) = 2;
    nDims(2) = 1;
    out = permute(in,nDims); 
end