classdef FeFunction < handle
    % NOTE
    % Go to P1Function.m for The Class Formerly Known As FeFunction
    % Eventually should extend Function/Field or something like that, to
    % account for other types of functions (eg. L2)
    
    properties (Constant, Access = public)
        fType = 'FE'
    end

    properties (Access = public)
       fValues        
    end

    properties (GetAccess = public, SetAccess = protected)
       ndimf
       order
       mesh       
    end
    
    properties (Access = protected)
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

        function n = computeL2norm(obj)
            l2Norm = L2Norm(obj.mesh);
            n = l2Norm.compute(obj);
        end
    end

    methods (Static, Access = public)
        function obj = create(type,fValues,mesh)
            s.order   = type;
            s.fValues = fValues;
            s.mesh    = mesh;
            obj       = LagrangianFunction(s);
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