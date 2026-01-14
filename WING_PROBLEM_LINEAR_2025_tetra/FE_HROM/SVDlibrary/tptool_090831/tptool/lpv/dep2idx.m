function depidx = dep2idx(dep)
%DEP2IDX Convert dependency data into indices (used internally)
%	depidx = DEP2IDX(dep)
%
%	dep    - dependency matrix for an LPV model (0:depends, 1:not depends)
%
%	See also SAMPLING_LPV

[O I P] = size(dep);
depidx = zeros([O I P]);
for i = 1:O
	for j = 1:I
		d = find(dep(i,j,:)>0);
		depidx(i,j,d) = 1:length(d);
	end
end
