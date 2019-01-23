function prova_ipopt

funcs = struct;
funcs.objective = @objective;
funcs.gradient = @gradient;
funcs.constraints = @constraints;
funcs.jacobian = @jacobian;
funcs.jacobianstructure = @jacobianstructure;

x0 = rand(4,1);
options.ipopt.print_level           = 5;
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.limited_memory_update_type = 'bfgs';
options.ub = ones(4,1);
options.lb = zeros(4,1);

[x, info] = ipopt(x0,funcs,options);

end

function f = objective (x)
    f = x(1)*x(4)*sum(x(1:3)) + x(3);
end

function g = gradient (x)
          g = [ x(1)*x(4) + x(4)*sum(x(1:3))
                x(1)*x(4)
                x(1)*x(4) + 1
                x(1)*sum(x(1:3)) ];
end

function c = constraints (x)
c = [ prod(x); sum(x.^2) ];
end

function J = jacobian (x)
    J = sparse([ prod(x)./x; 2*x ]);
end

function J = jacobianstructure() 
    J = sparse(ones(2,4));
end