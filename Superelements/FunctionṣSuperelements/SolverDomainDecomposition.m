function [GlobalU] = SolverDomainDecomposition(GlobalMatrices)
    LHS = GlobalMatrices.LHS;
    RHS = GlobalMatrices.RHS;

    s.type = "DIRECT";
    Sol = Solver.create(s);
    GlobalU = Sol.solve(LHS,RHS);
end