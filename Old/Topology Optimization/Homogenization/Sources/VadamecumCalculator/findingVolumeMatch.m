function findingVolumeMatch


vS = findVolume('SmoothRectangle');
vR = findVolume('Rectangle');










end


function volV = findVolume(micro)

d = load(['/home/alex/git-repos/SwanLab/Swan/Output/',micro,'/',micro,'.mat']);

var = d.d.variables;

mxV = d.d.domVariables.mxV;
myV = d.d.domVariables.myV;

for imx = 1:length(mxV)
    for imy = 1:length(myV)
        volV(imx,imy) = var{imx,imy}.volume;
    end
end



end