clc
clear
close all
% General call for testMicro --> General sampling
p.type = "SuperformulaInclusion";
% a = 0.01:0.01:0.499;
% b = 0.01:0.01:0.499;
% n = 0.5:0.1:35; 
p.a  = 0.25;
p.b  = 0.25;
p.m  = 5;
p.n1 = 2;
p.n2 = 6;
p.n3 = 6;
TestMicro(p);
% sample = 1;
% dataset = zeros(length(a)*length(b)*length(n),9);
% for i = 1:length(n)
%     for j = 1:length(a)
%         for k = 1:length(b)
%             gPar.a = a(j);
%             gPar.b = b(k);
%             gPar.n = n(i);
%             Ch = TestMicro(gPar).stateProblem.Chomog;
%             dataset(sample,:) = [gPar.a, gPar.b, gPar.n, Ch(1,1), Ch(1,2), Ch(1,3), Ch(2,2), Ch(2,3), Ch(3,3)];
%             sample = sample+1;
%         end
%     end
% end
% writematrix(dataset, "testDataset.csv", 'Delimiter', ",");

