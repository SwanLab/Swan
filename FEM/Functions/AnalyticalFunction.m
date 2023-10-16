classdef AnalyticalFunction < L2Function
    
    properties (Access = public)
        ndimf
    end
    
    properties (Access = private)
        fHandle
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = AnalyticalFunction(cParams)
            obj.init(cParams)
        end

        function fxV = evaluate(obj, xGLoc)
            xV = obj.mesh.computeXgauss(xGLoc);
            fxV = obj.fHandle(xV);
        end

        function plot(obj)
            p1D = obj.project('P1D');
            p1D.plot();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fHandle = cParams.fHandle;
            obj.ndimf   = cParams.ndimf;
            obj.mesh    = cParams.mesh;
        end
        
    end
    
end