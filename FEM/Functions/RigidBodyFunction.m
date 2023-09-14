classdef RigidBodyFunction < L2Function
    
    properties (Access = public)
        ndimf
    end
    
    properties (Access = private)
        fvalues
        refPoint
%         mesh

        fun
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = RigidBodyFunction(cParams)
            obj.init(cParams)
            obj.computeDisplacements()
        end

        function fxV = evaluate(obj, xGLoc)
            fxV = obj.fun.evaluate(xGLoc);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fvalues = cParams.fvalues;
            obj.refPoint = cParams.refPoint;
            obj.mesh    = cParams.mesh;
        end
       
        function computeDisplacements(obj)
            u=obj.fvalues(1);
            v=obj.fvalues(2);
            theta=obj.fvalues(3);
            x0 = obj.refPoint(1);
            y0 = obj.refPoint(2);
            s.fHandle=@(x) [u-(x(2,:,:)-y0)*theta ; v+(x(1,:,:)-x0)*theta];
            obj.ndimf   = obj.mesh.ndim;
            s.ndimf= obj.ndimf; 
            s.mesh=obj.mesh;
            obj.fun= AnalyticalFunction(s);
        end


    end
    
end