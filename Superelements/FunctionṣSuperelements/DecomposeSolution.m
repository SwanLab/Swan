function [Fields] = DecomposeSolution(DecomposeInputs)
    subK = DecomposeInputs.subMatrices.K;
    subC = DecomposeInputs.subMatrices.C;
    type = DecomposeInputs.type;
    super = DecomposeInputs.super;
    globalU = DecomposeInputs.globalU;
    
    if type == "Two"
        subdomains = size(subK,1);    
        boundaries = size(subC,1);
    elseif type == "Three"
        subD = DecomposeInputs.subMatrices.D;
        subdomains = size(subK,1);    
        boundaries = size(subC,1);
        interfaces = size(subD,1);
    end
    
    if super
         dofsBoundary = DecomposeInputs.dofsBoundary;
        [u,initialDof] = DecomposeSuperelementsDisplacement(globalU,subdomains,dofsBoundary);
        Fields.u = u;
    else
        [u,initialDof] = DecomposeDisplacement(globalU,subdomains,subK);
         Fields.u = u;
    end
    [lambda,initialDof] = DecomposeLagrange(globalU,boundaries,subC,initialDof,type);
    Fields.lambda = lambda;

    if type == "Three"
        [uI] = DecomposeInterfaceDisplacement(globalU,interfaces,subD,initialDof);
        Fields.uI = uI;
    end

   


end


