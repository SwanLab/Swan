classdef MinimumGradFieldWithVectorInL2 < handle

    properties (Access = private)
        LHS
        RHS
    end
    
    properties (Access = private)
        mesh
        fGauss
        dim
    end
    
    methods (Access = public)
        
        function obj = MinimumGradFieldWithVectorInL2(cParams)
            obj.init(cParams);
            obj.createDimension();
        end

        function u = solve(obj)
            obj.computeLHS();
            obj.computeRHS();
            u = obj.solveSystem();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh   = cParams.mesh;
            obj.fGauss = cParams.fGauss;
        end

        function createDimension(obj)
            q = Quadrature();
            q = q.set(obj.mesh.type);
            s.mesh = obj.mesh;
            s.pdim = '1D';
            s.ngaus = q.ngaus;
            d = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end
        
        function computeLHS(obj)
            K = obj.computeStiffnessMatrix();
            M = obj.computeMassMatrix();
            I = ones(size(K,1),1);
            %eta = 0.01;
            %obj.LHS = K + eta*M;
            obj.LHS = [K,I;I',0];
        end
        
        function K = computeStiffnessMatrix(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod        = obj.mesh.npnod;
            s.type         = 'StiffnessMatrix';
            s.dim          = obj.dim;
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
        end
        
        function M = computeMassMatrix(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod        = obj.mesh.npnod;
            s.type         = 'MassMatrix';
            s.dim          = obj.dim;
            s.quadType     = 'QUADRATIC';
            lhs = LHSintegrator.create(s);
            M = lhs.compute();
        end
        
        function computeRHS(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            s.fGauss    = obj.fGauss;
            s.xGauss    = q.posgp;
            s.mesh      = obj.mesh;
            s.type      = obj.mesh.type;
            s.quadOrder = q.order;
            rhs = RHSintegrator(s);
            rhsC = rhs.integrateWithShapeDerivative();
            rhsV = obj.assembleIntegrand(rhsC);
            obj.RHS = [rhsV;0];
        end
        
        function f = assembleIntegrand(obj,rhsCells)
            integrand = rhsCells;
            ndofs  = obj.mesh.npnod;
            connec = obj.mesh.connec;
            nnode  = size(connec,2);
            f = zeros(ndofs,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end
        
        function u = solveSystem(obj)
            s = Solver.create();
            u = s.solve(obj.LHS,obj.RHS);
            u = u(1:end-1);
        end
        
    end
    
end