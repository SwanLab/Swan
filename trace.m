function t = trace(A)

if size(A,1)~=size(A,2)
  error(message('MATLAB:trace:square'));
end

if ismatrix(A)
    t = full(sum(diag(A)));
else
    sum = 0;
    for i=1:size(A,1)
        diagElem = A(i,i,:,:);
        sum = sum + diagElem;
    end
    t = sum;
end
