function [diagonal,lower,names] = choleskyLayout(n)
[i,j]    = ndgrid(1:n,1:n);
diagonal = find(i == j);
lower    = find(i > j);
names    = ["logL" + string(i(diagonal)) + string(j(diagonal)); "L" + string(i(lower)) + string(j(lower))]';
end