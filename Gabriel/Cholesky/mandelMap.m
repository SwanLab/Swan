function A = mandelMap(names)
components = componentIndices(names);
n          = max(components(:));
nOut       = size(components,1);
A          = zeros(numel(toMandel(zeros(n,n,n,n))),nOut);
for k = 1:nOut
    C = zeros(n,n,n,n);
    c = components(k,:);
    for idx = [c; c([2 1 3 4]); c([1 2 4 3]); c([2 1 4 3]); c([3 4 1 2]); c([4 3 1 2]); c([3 4 2 1]); c([4 3 2 1])]'
        C(idx(1),idx(2),idx(3),idx(4)) = 1;
    end
    A(:,k) = reshape(toMandel(C),[],1);
end
end

function components = componentIndices(names)
names      = cellstr(names);
components = cell2mat(cellfun(@(name) name(2:end) - '0',names(:),'UniformOutput',false));
end