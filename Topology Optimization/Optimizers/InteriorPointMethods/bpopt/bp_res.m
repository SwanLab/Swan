% residual with equality constants and inequality slack variables
function [c] = bp_res(bp,x,s,bL,bU)

% get base residual
c = bp_res_stub(bp,x);

j = 0;
for i = 1:size(c,2),
    if (bU(i)==bL(i)),
        % equality constant
        c(i) = c(i) - bL(i);
    else
        % inequality slack
        j = j + 1;
        c(i) = c(i) - s(j);
    end
end
