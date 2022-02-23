classdef HarmonicProjection < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
       orientationAngle
       mesh        
    end
    
    properties (Access = private)
        dim
        massMatrix
        stiffnessMatrix
        rhsU
    end
    
    methods (Access = public)
        
        function obj = HarmonicProjection(cParams)
            obj.init(cParams);
            obj.createDimension();
            obj.computeMassMatrix();
            obj.computeStiffnessMatrix();
        end
        
        function project(obj)
            LHS = obj.computeLHS();
            rhs = obj.computeRHS(1);

        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.mesh             = cParams.mesh;
           obj.orientationAngle = cParams.orientationAngle;
        end
        
        function createDimension(obj)
            q = Quadrature();
            q = q.set(obj.mesh.type);
            s.mesh = obj.mesh;
            s.pdim = 'FILTER';
            s.ngaus = q.ngaus;
            d = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end        
        
        function computeMassMatrix(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod        = obj.mesh.npnod;
            s.type         = 'MassMatrix';
            s.dim          = obj.dim;
            lhs = LHSintegrator.create(s);
            M = lhs.compute();
            obj.massMatrix = M;
        end

       function K = computeStiffnessMatrix(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod        = obj.mesh.npnod;
            s.type         = 'StiffnessMatrix';
            s.dim          = obj.dim;
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
            obj.stiffnessMatrix = K;
       end           

       function LHS = computeLHS(obj)
           x = obj.mesh.coord(:,1);
           y = obj.mesh.coord(:,2);
           b = boundary(x,y,1);
           nInt = setdiff(1:obj.dim.npnod,b);
           K    = obj.stiffnessMatrix; 
           Kred = K(nInt,:);
           M    = obj.massMatrix;
           Z    = zeros(length(nInt),length(nInt));
           LHS  = [M,Kred';Kred,Z];
       end

        function computeRHS(obj,idim)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            s.mesh  = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod = obj.dim.npnod;
            s.dim   = obj.dim;
            s.type = 'SIMPLE';
            int = Integrator.create(s);
            rhsC = int.integrateFnodal(obj.orientationAngle(:,idim),q.order);
            obj.rhsU = obj.assembleIntegrand(rhsC);        
        end

    end
    
end