function [Ce] = fourier_law_elas(kxx,kxy,kyy)

Ce(1,1,:) = kxx;
Ce(1,2,:) = kxy;
Ce(2,1,:) = kxy;
Ce(2,2,:) = kyy;



end

