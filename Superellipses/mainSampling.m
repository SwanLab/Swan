clc
clear
close all
% General call for testMicro --> General sampling
gPar.type = "EllipseInclusion";
a = 0.2;
b = 0.2;
n = 0.2; 
% p.a  = 0.25;
% p.b  = 0.25;
% p.m  = 4;
% p.n1 = 1;
% p.n2 = 1;
% p.n3 = 1;
sample = 1;
dataset = zeros(length(a)*length(b)*length(n),9);
tic;
for i = 1:length(n)
    for j = 1:length(a)
        for k = 1:length(b)
            gPar.a = a(j);
            gPar.b = b(k);
            gPar.n = n(i);
            try
                Ch = TestMicro(gPar).stateProblem.Chomog;
                dataset(sample,:) = [gPar.a, gPar.b, gPar.n, Ch(1,1), Ch(1,2), Ch(1,3), Ch(2,2), Ch(2,3), Ch(3,3)];
                sample = sample+1;
                disp(sample)
            catch
                dataset(sample,:) = [];
            end
        end
    end
end
toc;
% writematrix(dataset, "finalSuperellipseDataset.csv", 'Delimiter', ",");

