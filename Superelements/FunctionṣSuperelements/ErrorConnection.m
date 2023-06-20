function [ErrorConnL2, ErrorConnLinf] = ErrorConnection(ErrorConnInputs)
    globalU = ErrorConnInputs.globalU;
    globalK = ErrorConnInputs.globalK;
    globalLHS = ErrorConnInputs.globalLHS;
    subC = ErrorConnInputs.subC;

    ErrorConnLinf = 0;
    ErrorConnL2 = 0;
    totalBoundaries = size(subC,1);
    initialDof = size(globalK,1);
    for i=1:totalBoundaries
        DofsBoundary = (1:1:+size(subC{i,1},2)) + initialDof;
        % Linf
        ErrorLinf = full(max(abs(globalLHS(DofsBoundary,:)*globalU)));
        if ErrorLinf > ErrorConnLinf
            ErrorConnLinf = ErrorLinf;
        end
        % L2
        ErrorL2 = full(sqrt(sum((globalLHS(DofsBoundary,:)*globalU).^2)));
        if ErrorL2 > ErrorConnL2
            ErrorConnL2 = ErrorL2;
        end
        initialDof = max(DofsBoundary);
    end
end