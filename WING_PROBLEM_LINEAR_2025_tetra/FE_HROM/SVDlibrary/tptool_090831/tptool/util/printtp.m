function printtp(S)
%PRINTTP Print LTI models of a TP model
%	PRINTTP(S)
%	
%	S - core tensor of a TP model

n = size(S);
LTI = shiftdim(S, length(n)-2)
