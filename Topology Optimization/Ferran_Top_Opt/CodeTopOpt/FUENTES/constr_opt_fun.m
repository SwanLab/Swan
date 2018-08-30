function [constr,ceq,gradient,geq] = constr_opt_fun (x,constr_cell) 

nconstr = size(constr_cell,1);
nvar = size(x,1);

constr = zeros(nconstr,1);
gradient = zeros(nvar,nconstr);


for i = 1:nconstr
    [constr(i),gradient(:,i)] = constr_cell{i}(x);
end

ceq = [];
geq = [];

end