function u = vcycle(u, b, data, vdown, vup, level, bc, numero, mesh)

    if level == 1
        u = data(level).A\b;
    else 
        %[u,~] = gauss_seidel(data(level).A, b, u,vdown);
        [u, res] = conjugateGradient_Solver(data(level).A,b,u);
        malla = 'fine';
        %plotSolution(u, mesh{1,2}, bc{1,2}, malla, numero)
        %plotRes(res,mesh{1,2},bc{1,2},malla,numero)
        numero = numero + 1;
        r = b - data(level).A*u;
        r = bc{level}.reducedToFullVector(r);
        r = reshape(r,2,[])';
        Rr = data(level - 1).R * r;
        Rr = reshape(Rr,1,[])';
        Rr = bc{level - 1}.fullToReducedVector(Rr);
        e = vcycle(0*Rr, Rr, data, vdown, vup, level - 1);
        e = bc{level - 1}.reducedToFullVector(e);
        e = reshape(e,2,[])';
        e = data(level - 1).T * e; 
        e = reshape(e,1,[])';
        e = bc{level}.fullToReducedVector(e);
        u = u + e;
        %[u,~] = gauss_seidel(data(level).A, b, u,vup);
        [u, res] = conjugateGradient_Solver(data(level).A,b,u);
        malla = 'coarse';
        %plotSolution(u, mesh{1,1}, bc{1,1}, malla, numero)
        %plotRes(res,mesh{1,1},bc{1,1},malla,numero)
        numero = numero + 1;
    end

end