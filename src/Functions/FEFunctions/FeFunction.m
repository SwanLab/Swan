classdef FeFunction < BaseFunction
    
    properties (Constant, Access = public)
        fType = 'FE'
    end

    properties (GetAccess = public, SetAccess = protected)
        fValues
        order      
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