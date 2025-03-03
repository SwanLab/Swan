function findingSampleWithMaxVolume

d = load('/home/alex/git-repos/SwanLab/Swan/Output/Rectangle/Rectangle.mat');

var = d.d.variables;

mxV = d.d.domVariables.mxV;
myV = d.d.domVariables.myV;

for imx = 1:length(mxV)
    for imy = 1:length(myV)
        volV(imx,imy) = var{imx,imy}.volume;
    end
end

N = 10000;
mx = linspace(min(mxV),max(mxV),N);
my = linspace(min(myV),max(myV),N);
vol = interp2(mxV,myV,volV,mx,my);

for imx = 1:length(mx)
    v(imx) = vol(imx);
end

[vmin,imax] = min(abs(v-0.0939));


vol(imax)
mx(imax)
my(imax)
number = length(mxV)*(imax-1) + imax

end
