function [fun, funNames] = PrintTensor(tensor,mesh,filename)
    basename = inputname(1);
    dof = 0;
    xV = Quadrature.create(mesh,1).posgp;
    tens_eval = squeeze(tensor.evaluate(xV));
    fun = {};
    funNames = {};
    for iComp = 1:size(tens_eval,1)
        for jComp = 1:size(tens_eval,2)
            dof = dof + 1;
            a.fValues = tens_eval(iComp, jComp, :);
            a.mesh   = mesh;
            a.order  = 'P0';
            p0 = LagrangianFunction(a);
            fun{dof} = p0;
            funNames{dof} = strcat(basename,string(iComp),string(jComp));
        end
    end
    b.mesh     = mesh;
    b.filename = filename;
    b.fun      = fun;
    b.funNames = funNames;
    b.type     = 'Paraview';
    pst = FunctionPrinter.create(b);
    pst.print();


end