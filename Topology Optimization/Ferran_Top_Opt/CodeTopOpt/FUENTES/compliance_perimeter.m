function [ f,df ] = compliance_perimeter( x,compliance,perimeter,lambda )

[c,dc] = compliance(x);
[p,dp] = perimeter(x);

f = c + lambda*p;
df = dc + lambda*dp;

end

