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
                dimC       = ones(1,ndims(constant));
                s.constant = constant;
                s.ndimf    = size(constant);
                s.mesh     = mesh;
                s.fHandle = @(xV) squeezeParticular(repmat(constant,[dimC,size(xV,2),mesh.nelem]),2);
                obj = ConstantFunction(s);
            end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fHandle = cParams.fHandle;
            obj.constant = cParams.constant;
            obj.mesh = cParams.mesh;
            obj.ndimf = cParams.ndimf;
            obj.ndimfTotal = prod(obj.ndimf);
        end

    end

    methods (Access = protected)

        function fxV = evaluateNew(obj, xV)
            fxV = obj.fHandle(xV);
        end

    end

end