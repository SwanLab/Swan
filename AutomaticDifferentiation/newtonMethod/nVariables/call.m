clear;
nodes = 3;
u = 1 + zeros(1,nodes);
[plotU,plotVal,plotGrad,plotGrad2,iterations] = newtonMethodFindMinNVariables(u);