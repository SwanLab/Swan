    
aEval = [1 2;3 4];
aEval = repmat(aEval,1,1,2,2,2);
bEval = [0 5; 6 7];
bEval = repmat(bEval,1,1,2,2,2);

sizeA = size(aEval);
sizeB = size(bEval);
if sizeA(3:end) == sizeB(3:end)
    res = zeros([sizeA(1)*sizeB(1),sizeA(2)*sizeB(2),sizeA(3:end)]);
    for i=1:sizeA(1)
        for j=1:sizeA(2)
            for k=1:sizeB(1)
                for l=1:sizeB(2)
                    indexA = [i,j,repmat({':'}, 1, ndims(aEval)-2)];
                    Aij = squeeze(aEval(indexA{:}));
                    indexB = [k,l,repmat({':'}, 1, ndims(bEval)-2)];
                    Bkl = squeeze(bEval(indexB{:}));

                    idx1 = sizeB(1)*(i-1) + k;
                    idx2 = sizeB(2)*(j-1) + l;
                    indexRes = [idx1,idx2,repmat({':'}, 1, ndims(res)-2)];
                    res(indexRes{:}) = Aij.*Bkl;
                end
            end
        end
    end
    fVR = res;
else
    error('Array dimensions are not the same')
end