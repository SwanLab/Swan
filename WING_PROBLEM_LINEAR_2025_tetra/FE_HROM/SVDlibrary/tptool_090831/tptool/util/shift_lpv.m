function S = shift_lpv(D, v)
%SHIFT_LPV Shift a sampled LPV model by a vertex system
%	S = SHIFT_LPV(D, v)
%
%	D - sampled LPV system
%	v - shift matrix
%
%	S - shifted LPV system data


S = D;
[Sy Sx] = size(D);

for i = 1:Sy
	for j = 1:Sx
		S{i,j} = S{i,j} + v(i,j);
	end
end
