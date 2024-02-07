function [u,numero,residuFine,iterRes] = vcycle(u, b, data, vdown, vup, level, bc, numero, mesh,residuFine,iterRes)

    if level == 1
        %u = data(level).A\b;
        malla = 'coarse';
        [u, res,residuFine,iterRes] = conjugateGradient_Solver(data(level).A,b,u,malla,residuFine,iterRes);
        %plotSolution(u, mesh{1,1}, bc{1,1}, malla, numero)
        %plotRes(res,mesh{1,1},bc{1,1},malla,numero)
        numero = numero + 1;
    else 
        %[u,~] = gauss_seidel(data(level).A, b, u,vdown);
        malla = 'fine';
        [u, res,residuFine,iterRes] = conjugateGradient_Solver(data(level).A,b,u,malla,residuFine,iterRes);
        %plotSolution(u, mesh{1,2}, bc{1,2}, malla, numero)
        %plotRes(res,mesh{1,2},bc{1,2},malla,numero)
        numero = numero + 1;
        r = b - data(level).A*u;
        r = bc{level}.reducedToFullVector(r);
        r = reshape(r,2,[])';
        Rr = data(level - 1).R * r;
        Rr = reshape(Rr,1,[])';
        Rr = bc{level - 1}.fullToReducedVector(Rr);
        [e,numero,residuFine,iterRes] = vcycle(0*Rr, Rr, data, vdown, vup, level - 1, bc, numero, mesh,residuFine,iterRes);
        %[e,numero] = vcycle(Rr, data(level - 1).b, data, vdown, vup, level - 1, bc, numero, mesh);
        e = bc{level - 1}.reducedToFullVector(e);
        e = reshape(e,2,[])';
        e = data(level - 1).T * e; 
        e = reshape(e,1,[])';
        e = bc{level}.fullToReducedVector(e);
        u = u + e;
        %[u,~] = gauss_seidel(data(level).A, b, u,vup);
        malla = 'fine';
        [u, res,residuFine,iterRes] = conjugateGradient_Solver(data(level).A,b,u,malla,residuFine,iterRes);
        %plotSolution(u, mesh{1,2}, bc{1,2}, malla, numero)
        %plotRes(res,mesh{1,2},bc{1,2},malla,numero)
        numero = numero + 1;
     end

end