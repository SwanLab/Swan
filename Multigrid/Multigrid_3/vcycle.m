function u = vcycle(u, b, data, vdown, vup, level)
if level == 1
    u = data(level).A\b;
else 
    [u,~] = gauss_seidel(data(level).A, b, u,vdown);
    r = b - data(level).A*u;
    Rr = data(level - 1).R * r;
    e = data(level - 1).T * vcycle(0*Rr, Rr, data, vdown, vup, level - 1); 
    u = u + e;
    [u,~] = gauss_seidel(data(level).A, b, u,vup);
end

end