classdef HarmonicProjection < handle
    
    properties (Access = private)
        dim
        massMatrix
        stiffnessMatrix
        reducedStiffnessMatrix
        LHS
        solver
    end    

    properties (Access = private)
       mesh   
       boundaryMesh
    end
    
    methods (Access = public)
        
        function obj = HarmonicProjection(cParams)
            obj.init(cParams);
            obj.createDimension();
            obj.computeMassMatrix();
            obj.computeStiffnessMatrix();
            obj.computeReducedStiffnessMatrix();
            obj.computeLHS();
            obj.createSolver()
        end
        
        function [vH,errF] = project(obj,v)
            lhs = obj.LHS;
            rhs = obj.computeRHS(v);
            vH  = obj.solver.solve(lhs,rhs);
            vH  = vH(1:obj.dim.npnod,1);

            Kred = obj.reducedStiffnessMatrix;            
            M    = obj.massMatrix;
            grad = Kred*vH;
            errC = norm(grad);
            errF = (v-vH)'*M*(v-vH)/(v'*M*v);
        end

        function [vH,errF] = projectByDual(obj,v)
            Kred = obj.reducedStiffnessMatrix;            
            M    = obj.massMatrix;            
            lhs = Kred*(M\Kred');
            rhs = Kred*v;            
            lambda  = obj.solver.solve(lhs,rhs);
            vH = v - M\(Kred'*lambda);
            errC = norm(Kred*vH);
            errF = (v-vH)'*M*(v-vH)/(v'*M*v);
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.mesh             = cParams.mesh;
           obj.boundaryMesh     = cParams.boundaryMesh;
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
            %M = diag(sum(M));
            %M = eye(size(M,1));
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

       function computeReducedStiffnessMatrix(obj)
           b    = obj.boundaryMesh;
           nInt = setdiff(1:obj.dim.npnod,b);
           K    = obj.stiffnessMatrix; 
           Kred = K(nInt,:);
           obj.reducedStiffnessMatrix = Kred;
       end

       function createSolver(obj)
            s = Solver.create();
            obj.solver = s;
       end

       function  computeLHS(obj)
           Kred = obj.reducedStiffnessMatrix;
           M    = obj.massMatrix;
           Z    = obj.computeZeroFunction();
           lhs  = [M,Kred';Kred,Z];
           obj.LHS = lhs;
       end

       function Z = computeZeroFunction(obj)
           b    = obj.boundaryMesh;
           nInt = setdiff(1:obj.dim.npnod,b);           
           Z    = zeros(length(nInt),length(nInt));           
       end

        function rhs = computeRHS(obj,v)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            s.mesh  = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod = obj.dim.npnod;
            s.dim   = obj.dim;
            s.type = 'SIMPLE';
            int = Integrator.create(s);
            rhs = int.integrateFnodal(v,q.order);
            b = obj.boundaryMesh;
            nInt = setdiff(1:obj.dim.npnod,b);
            Z   = zeros(length(nInt),1);
            rhs = [rhs;Z];
        end

    end
    
end