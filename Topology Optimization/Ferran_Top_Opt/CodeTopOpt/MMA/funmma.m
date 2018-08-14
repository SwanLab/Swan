function [f,df,c,dc] = funmma(x,fobj,fconstr,constraint_type)

[f,df] = fobj(x);
[c,~,dc,~] = fconstr(x);
dc = dc';

%% Check constraint case
if strcmp(constraint_type,'equality')
    c = [c;-c];
    dc = [dc;-dc];
end

%% Re-scale constraints
% In many applications, the constraints are on the form yi(x) <= ymaxi 
% The user should then preferably scale the constraints in such a way that 1 <= ymaxi <= 100 for each i 
% (and not ymaxi = 10^10 for example).
kconstr = 1;
cconstr = 0;
c = kconstr*c;
c(c > 0) = c(c > 0) + cconstr;
c(c < 0) = c(c < 0) - cconstr;
% dc = kconstr*dc;

%% Re-scale objective function
% The objective function f(x) should preferably be scaled such that
% 1 <= f0(x) <= 100 for reasonable values on the variables.
kfun = 1;
cfun = 0;
f = kfun*f + cfun;
df = kfun*df;

end