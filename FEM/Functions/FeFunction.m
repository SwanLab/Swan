classdef FeFunction < handle
    % NOTE
    % Go to P1Function.m for The Class Formerly Known As FeFunction
    % Eventually should extend Function/Field or something like that, to
    % account for other types of functions (eg. L2)
    
    properties (Constant, Access = public)
        fType = 'FE'
    end

    properties (Access = public)
       ndimf
       order
       fValues
    end
    
    properties (Access = protected)
        mesh
    end

    properties (Access = private)
        
    end
    
    methods (Access = public)
        function fun = project(obj,target)
            s.mesh          = obj.mesh;
            s.projectorType = target;
            proj = Projector.create(s);
            fun = proj.project(obj);
        end

        function totVal = computeScalarProduct(obj,f,order)
            q.mesh     = obj.mesh;
            q.quadType = order;
            q.type     = 'ScalarProduct';
            int        = Integrator.create(q);
            totVal     = int.compute(obj,f);
        end
    end

    methods (Static, Access = public)
        function obj = create(cParams)
            fun = FunctionFactory();
            obj = fun.create(cParams);
        end

        function obj = createEmpty(cParams)
            feFunType = cParams.feFunType;
            mesh      = cParams.mesh;
            ndimf     = int2str(cParams.ndimf);
            specs     = ['.create(mesh,',ndimf,')'];
            obj       = eval([feFunType,specs]);
        end
    end

end