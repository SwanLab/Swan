function [u, res_record]  = mgsolve(data, vdown, vup,tol)
nlevels = size(data,2)-1;
res = data(nlevels).b;
i = 2;
u = 0*data(nlevels).b;
res_record = norm(data(nlevels).b,inf);

    while norm(data(nlevels).b - data(nlevels).A * u, inf) >= tol 
       u = vcycle(u, data(nlevels).b, data, vdown, vup, nlevels);
       res_record(i) = norm(data(nlevels).b - data(nlevels).A * u,inf);
       i = i + 1;
    end
    
end