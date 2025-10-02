function dom = VonMisesStress(stress)
    s.operation = @(xV) evaluate(stress,xV);
    s.mesh      = stress.mesh;
    dom         = DomainFunction(s);
end

function fEval = evaluate(stress, xV)
    m = stress.mesh;
    for i=1:m.ndim
        for j=1:m.ndim
        base   = createBaseFunction(m,i,j);
        stressMatrix = DDP(stress,base);
        sigma{i,j} = stressMatrix.evaluate(xV);
        end    
    end

    switch m.ndim
        case 2
            fEval = sqrt( sigma{1,1}.^2 + sigma{2,2}.^2 ...
                - sigma{1,1}.*sigma{2,2} ...
                + 3*sigma{1,2}.^2 );
        case 3
          fEval = sqrt(0.5*((sigma{1,1}-sigma{2,2}).^2 + (sigma{2,2}-sigma{3,3}).^2 ...
                + (sigma{3,3}-sigma{1,1}).^2 )+3*(sigma{1,2}.^2 + sigma{2,3}.^2 + sigma{1,3}.^2));
    end
end

function base = createBaseFunction(m,i,j)
    constant = zeros(2,2);
    constant(i,j) = 1;
    base = ConstantFunction.create(constant,m);
end