function [ c,dc ] = constraint_CCstar( x,compliance,epsilon )

[c,dc] = compliance(x);

c = c - epsilon;

end