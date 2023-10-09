classdef RigidBodyFunction < L2Function
    
    properties (Access = public)
        ndimf
    end
    
    properties (Access = private)
        fvalues
        refPoint
%         mesh

        fun
        
        horizontalTranslationBase
        verticalTranslationBase
        rotationalBase
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = RigidBodyFunction(cParams)
            obj.init(cParams)
            obj.computeHorizontalTranslationBase();
            obj.computeVerticalTranslationBase();
            obj.computeRotationBase();
        end

        function fxV = evaluate(obj, xGLoc)
            phiU = obj.horizontalTranslationBase.evaluate(xGLoc);
            phiV = obj.verticalTranslationBase.evaluate(xGLoc);
            phiT = obj.rotationalBase.evaluate(xGLoc);
            u     = obj.fvalues(1);
            v     = obj.fvalues(2);
            theta = obj.fvalues(3);             
            fxV = u*phiU + v*phiV + theta*phiT;
        end
        
        function basis = computeBasisFunction(obj,xGloc)
            basis{1} = obj.horizontalTranslationBase.evaluate(xGloc);
            basis{2} = obj.verticalTranslationBase.evaluate(xGloc);
            basis{3} = obj.rotationalBase.evaluate(xGloc);
        end

    end

    methods (Access = public, Static)

        function RB = create(mesh,refPoint)
            ndimf=mesh.ndim;
            s.fvalues  = zeros(mesh.nnodes, ndimf);
            s.mesh     = mesh;
            s.refPoint = refPoint;
            RB = RigidBodyFunction(s);
        end

        function fS = times(f1,f2)
            fS = f1.fValues.*f2.fValues;
            s.fValues = fS;
            s.mesh    = f1.mesh;
            fS = RigidBodyFunction(s);
        end

        function fS = sum(f1,f2)
            fS = f1.fValues+f2.fValues;
            s.fValues = fS;
            s.mesh    = f1.mesh;
            fS = RigidBodyFunction(s);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fvalues  = cParams.fvalues;
            obj.refPoint = cParams.refPoint;
            obj.mesh     = cParams.mesh;
        end

        function computeHorizontalTranslationBase(obj)
            s.fHandle = @(x) [ones(size(x(1,:,:)));zeros(size(x(1,:,:)))];
            obj.ndimf = obj.mesh.ndim;
            s.ndimf   = obj.ndimf; 
            s.mesh    = obj.mesh;
            obj.horizontalTranslationBase = AnalyticalFunction(s);
        end    

        function computeVerticalTranslationBase(obj)
            s.fHandle = @(x) [zeros(size(x(1,:,:)));ones(size(x(2,:,:)))];
            obj.ndimf = obj.mesh.ndim;
            s.ndimf   = obj.ndimf; 
            s.mesh    = obj.mesh;
            obj.verticalTranslationBase = AnalyticalFunction(s);
        end              

        function computeRotationBase(obj)
            x0 = obj.refPoint(1);
            y0 = obj.refPoint(2);            
            s.fHandle = @(x) [-(x(2,:,:)-y0);x(1,:,:)-x0];
            obj.ndimf = obj.mesh.ndim;
            s.ndimf   = obj.ndimf; 
            s.mesh    = obj.mesh;
            obj.rotationalBase = AnalyticalFunction(s);
        end                
     

    end
    
end