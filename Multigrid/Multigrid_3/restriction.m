function R = restriction(T)
%computes Restriction matrices for a piecewise linear interpolant\
T = full(T);
R = T.';

for i = 1:size(R,1)
    if sum(R(i,:)) ~= 1
        R(i,1:end) = R(i,1:end)./sum(R(i,:));           
    end
end

end

