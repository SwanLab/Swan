clear;
clc;
close all;

filename = 'anisoCantilever';
a.fileName = filename;
gid = FemDataContainer(a);
mesh = gid.mesh;

s.fHandle = @(x) 1-heaviside((x(1,:,:)-1).^2+(x(2,:,:)-0.5).^2-0.3.^2);
s.ndimf   = 1;
s.mesh    = mesh;
fun{1}    = AnalyticalFunction(s);

s.mesh  = mesh;
s.alpha = 4;
s.beta  = 1;
s.theta = 90;
filter  = NonLinearFilterSegment(s);
filter.updateEpsilon(1);
[fun{2},err] = filter.compute(fun{1},2);

betaVec = 0.99:-0.01:0;
for i = 1:length(betaVec)
    filter.updateBeta(betaVec(i));
    [fun{end+1},newErr] = filter.compute(fun{1},2);
    err = [err;newErr];
end

figure
plot(err)
fun{end}.plot();