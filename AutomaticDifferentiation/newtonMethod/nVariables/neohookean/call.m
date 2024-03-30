clear;
nodes = 3;
u0 = 0.1 + zeros(1,nodes);
[plotU,plotVal,plotGrad,plotGrad2,iterations] = newtonMethodFindMinNVariablesNeohookean(u0);