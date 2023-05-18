clear;
[plotU,plotVal,plotGrad,plotGrad2,iterations] = newtonMethodFindMin([1 1 1]);
disp(min(plotU));