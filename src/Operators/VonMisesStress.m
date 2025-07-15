function dom = VonMisesStress(stress)
    s.operation = @(xV) evaluate(stress,xV);
    s.mesh      = stress.mesh;
    dom         = DomainFunction(s);
end

function fEval = evaluate(stress, xV)
    m = stress.mesh;
    for i=1:3
        for j=1:3 
        base   = createBaseFunction(m,i,j);
        stressMatrix = DDP(stress,base);
        sigma{i,j} = stressMatrix.evaluate(xV);
        end    
    end
    fEval = sqrt(0.5*((sigma{1,1}-sigma{2,2}).^2 + (sigma{2,2}-sigma{3,3}).^2 ...
                + (sigma{3,3}-sigma{1,1}).^2 )+3*(sigma{1,2}.^2 + sigma{2,3}.^2 + sigma{1,3}.^2));
end

function base = createBaseFunction(m,i,j)
    constant = zeros(3,3);
    constant(i,j) = 1;
    base = ConstantFunction.create(constant,m);
end