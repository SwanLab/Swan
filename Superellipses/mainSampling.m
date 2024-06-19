clc
clear
close all
% General call for testMicro --> General sampling
gPar.type = "EllipseInclusion";
gPar.a  = 0.37;
gPar.b  = 0.37;
gPar.n  = 2;
sample = 1;
dataset = [];
tic;
try
    Ch = TestMicro(gPar).stateProblem.Chomog;
    dataset(sample,:) = [gPar.a,gPar.b,gPar.n, Ch(1,1), Ch(1,2), Ch(1,3), Ch(2,2), Ch(2,3), Ch(3,3)];
    sample = sample+1;
    disp(sample)
catch
    dataset(sample,:) = [];
end
toc;
% writematrix(dataset, "finalEllipseDataset.csv", 'Delimiter', ",");

