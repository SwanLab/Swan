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

        function r = times(obj,b)
            s.operation = @(xV) obj.evaluate(xV);
            f           = DomainFunction(s);
            r           = f.*b;
        end
        
        function ord = getOrderNum(obj)
            ord = 2;
        end

    end

    methods (Access = public, Static)
        
        function obj = create(fHandle,ndimf,mesh)
            s.fHandle = fHandle;
            s.ndimf   = ndimf;
            s.mesh    = mesh;
            obj = AnalyticalFunction(s);
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
