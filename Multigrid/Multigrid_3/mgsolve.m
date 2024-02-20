function [u, res_record]  = mgsolve(data, vdown, vup, tol, bc, mesh)
%nlevels = size(data,2)-1;
nlevels = size(data,2);
%res = data(nlevels).b;
i = 1;
u = 0*data(nlevels).b;
res_record = norm(data(nlevels).b,inf);
numero = 1;
residuFine = 0;
iterRes = 1;

    while norm(data(nlevels).b - data(nlevels).A * u, inf) >= tol 
       [u,numero,residuFine,iterRes] = vcycle(u, data(nlevels).b, data, vdown, vup, nlevels, bc, numero, mesh, residuFine,iterRes);
       res_record(i) = norm(data(nlevels).b - data(nlevels).A * u,inf);
       i = i + 1;
    end

    res = residuFine;
    iter = linspace(1, length(residuFine), length(residuFine));
    plot(iter,log(residuFine))
    
end