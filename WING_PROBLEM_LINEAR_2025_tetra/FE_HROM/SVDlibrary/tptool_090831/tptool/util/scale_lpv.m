function S = scale_lpv(D, v)
%SCALE_LPV Scale each element of an LPV model
%	S = SCALE_LPV(D, v)
%
%	D - sampled LPV model
%	v - scale matrix for each element
%
%	S - scaled LPV model

S = D;
[Sy Sx] = size(D);
for i = 1:Sy
	for j = 1:Sx
		D{i,j} = D{i,j} * v(i,j);
	end
end
