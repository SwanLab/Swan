classdef FeFunction < BaseFunction
    
    properties (Constant, Access = public)
        fType = 'FE'
    end

    properties (Access = public)
       fValues        
    end

    properties (GetAccess = public, SetAccess = protected)
       order      
    end
    
    properties (Access = protected)
    end

    properties (Access = private)
        
    end
    
    methods (Access = public)


        

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
            ndimf     = int2str(cParams.ndimf);
            specs     = ['.create(mesh,',ndimf,')'];
            obj       = eval([feFunType,specs]);
        end
    end

end