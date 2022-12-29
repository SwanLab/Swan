classdef MinimumDiscGradFieldWithVectorInL2 < handle

    properties (Access = private)
        LHS
        RHS
    end
    
    properties (Access = private)
        mesh
        rhsType
        meshCont
        fGauss
        interp
        field
        interpolator
    end
    
    methods (Access = public)
        
        function obj = MinimumDiscGradFieldWithVectorInL2(cParams)
            obj.init(cParams);
            obj.createField();
        end

        function u = solve(obj)
            obj.computeLHS();
            obj.computeRHS();
            uC = obj.solveSystem();
            In = obj.interpolator; 
            u  = In*uC; 
            u = reshape(u,obj.mesh.nnodeElem,[])'; % Eh          
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.meshCont     = cParams.mesh;
            obj.rhsType      = cParams.rhsType;
            obj.mesh         = obj.meshCont.createDiscontinuousMesh();            
            obj.fGauss       = cParams.fGauss;
            obj.interpolator = cParams.interpolator;
        end

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = obj.mesh.interpolation.order;
            obj.field = Field(s);
        end

%         
        function computeLHS(obj)
            K = obj.computeStiffnessMatrix();
            In = obj.interpolator;
            K = In'*K*In;
           % M = obj.computeMassMatrix();
            %I = ones(size(K,1),1);
            %eta = 0.01;
            %obj.LHS = K + eta*M;
            %obj.LHS = [K,I;I',0];
            obj.LHS = K;
        end
        
        function K = computeStiffnessMatrix(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.type         = 'StiffnessMatrix';
            s.field        = obj.field;

            lhs = LHSintegrator.create(s);
            K = lhs.compute();
        end
        
        function computeRHS(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            s.fType     = 'Gauss';
            s.fGauss    = obj.fGauss;
            s.xGauss    = q.posgp;
            s.mesh      = obj.mesh;
            s.type      = obj.mesh.type;
            s.quadOrder = q.order;
            s.npnod     = obj.field.dim.ndofs;
            s.type      = obj.rhsType;
            s.globalConnec = obj.mesh.connec;
            rhs  = RHSintegrator.create(s);
            rhsV = rhs.compute();        
            In = obj.interpolator;
            rhsV = In'*rhsV;
            %obj.RHS = [rhsV;0];
            obj.RHS = rhsV;
        end
        
        function f = assembleIntegrand(obj,rhsCells)
            integrand = rhsCells;
            ndofs  = obj.mesh.nnodes;
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
            a.type = 'DIRECT';
            s = Solver.create(a);
            u = s.solve(obj.LHS,obj.RHS);
            u = u(1:end);
        end
        
    end
    
end