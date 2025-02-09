function gaussian = gaussianFunction
s = 0.1;
gaussian = @(txi,psi) 1/sqrt(2*pi*s^2)*exp(-(txi - psi).^2/(2*s^2));
end