classdef ConstantFunction < BaseFunction
    
    properties (GetAccess = public, SetAccess = private)
        constant
    end

    properties
        fHandle
    end
    
    
    methods (Access = public)
        
        function obj = ConstantFunction(cParams)
            obj.init(cParams)
        end
        
    end
    
    methods (Access = public, Static)
            
            function obj = create(constant, mesh)
                s.constant = constant;
                s.ndimf = length(constant);
                s.mesh = mesh;
                s.fHandle = @(xV) constant.*ones([length(constant),size(xV,2),mesh.nelem]); 
                obj = ConstantFunction(s);
            end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fHandle = cParams.fHandle;
            obj.constant = cParams.constant;
            obj.mesh = cParams.mesh;
            obj.ndimf = cParams.ndimf;
        end

    end

    methods (Access = protected)

        function fxV = evaluateNew(obj, xV)
            fxV = obj.fHandle(xV);
        end

    end

end