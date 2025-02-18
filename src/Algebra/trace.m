function t = trace(A)

if size(A,1)~=size(A,2)
  error(message('MATLAB:trace:square'));
end

if ismatrix(A)
    t = full(sum(diag(A)));
else
    t = 0;
    for i=1:size(A,1)
        diagElem = squeezeParticular(A(i,i,:,:),1);
        t = t + diagElem;
    end
end