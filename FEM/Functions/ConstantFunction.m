classdef ConstantFunction < AnalyticalFunction
    
    properties (GetAccess = public, SetAccess = private)
        constant
    end
    
    
    methods (Access = public)
        
        function obj = ConstantFunction(cParams)
            obj@AnalyticalFunction(cParams)
            obj.constant = cParams.constant;     
        end
        
    end
    
    methods (Access = public, Static)
            
            function obj = create(constant, mesh)
                s.constant = constant;
                s.ndimf = length(constant);
                s.mesh = mesh;
                s.fHandle = @(xV) constant*ones([1,size(xV,[2,3])]); 
                obj = ConstantFunction(s);
            end
    end
    
end