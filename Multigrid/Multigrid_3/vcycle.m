function [u,numero] = vcycle(u, b, data, vdown, vup, level, bc, numero, mesh, meshType)

    if level == 1
        %u = data(level).A\b;
        meshType = 'coarse';
        maxIter = 1000;
        %plotSolution(u, mesh{1,1}, bc{1,1}, malla, numero)
        [u, res] = conjugateGradient_Solver(data(level).A,b,u,meshType,maxIter);
        %plotSolution(u, mesh{1,1}, bc{1,1}, meshType, numero)
        %plotRes(res,mesh{1,1},bc{1,1},meshType,numero)
        numero = numero + 1;
    else 
        %[u,~] = gauss_seidel(data(level).A, b, u,vdown);
        meshType = 'fine';
        maxIter = 20;
        [u, res] = conjugateGradient_Solver(data(level).A,b,u,meshType,maxIter);
        %plotSolution(u, mesh{1,2}, bc{1,2}, meshType, numero)
        %plotRes(res,mesh{1,2},bc{1,2},meshType,numero)
        numero = numero + 1;
        r = b - data(level).A*u;
        [Rr] = interpolate(r,bc,data,level);
        [ur] = interpolate(u,bc,data,level);
        [er,numero] = vcycle(0*ur, Rr, data, vdown, vup, level - 1, bc, numero, mesh, meshType);
        e = restriction(er,bc,data,level);
        u = u + e;
        %[u,~] = gauss_seidel(data(level).A, b, u,vup);
        meshType = 'fine';
        [u, res] = conjugateGradient_Solver(data(level).A,b,u,meshType,maxIter);
        %plotSolution(u, mesh{1,2}, bc{1,2}, meshType, numero)
        %plotRes(res,mesh{1,2},bc{1,2},meshType,numero)
        numero = numero + 1;
     end

end

function fFine = restriction(fCoarse,bc,data,level)
        fFine = bc{level - 1}.reducedToFullVector(fCoarse);
        fFine = reshape(fFine,2,[])';
        fFine = data(level - 1).T * fFine; 
        fFine = reshape(fFine',[],1);
        fFine = bc{level}.fullToReducedVector(fFine);
end

function [fCoarse] = interpolate(fFine,bc,data,level)
        
        fFine = bc{level}.reducedToFullVector(fFine);
        fFine = reshape(fFine,2,[])';
        fCoarse = data(level - 1).R * fFine;
        fCoarse = reshape(fCoarse',[],1);
        fCoarse = bc{level - 1}.fullToReducedVector(fCoarse);
end
