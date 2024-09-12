function [ nodfun ] = mixedpdeprtni( p,t,elfun )
%mixedpdeprtni Element to nodal values
n = size(elfun,1);

for i=1:n
    [~,unitM,F] = assema(p,t,0,1,elfun(i,:));
    nodfun(:,i) = unitM\F;
end
end

