clear all
n = 50;
axV = linspace(0.5,0.99,n);
ayV = linspace(0.5,0.99,n);

for ix = 1:n
    for iy = 1:n
        ax = axV(ix);
        ay = ayV(iy);
        r(iy,ix) = (1-ay)/(1-ax)*(ax/ay);
    end
end
surf(axV,ayV,r)
contour(r,500)