classdef MinimumDiscGradFieldWithVectorInL2 < handle

    properties (Access = private)
        LHS
        RHS
    end
    
    properties (Access = private)
        mesh
        meshCont
        fGauss
        fValues
        field
        fieldDisc
        symCond
    end
    
    methods (Access = public)
        
        function obj = MinimumDiscGradFieldWithVectorInL2(cParams)
            obj.init(cParams);
            obj.createField();
        end

        function u = solve(obj)
            obj.createDiscontinousField();
            obj.createSymmetricMapCond()
            obj.computeLHS();
            obj.computeRHS();
            u = obj.solveSystem();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.meshCont = cParams.mesh;
            obj.mesh     = obj.meshCont.createDiscontinousMesh();
            obj.fGauss   = cParams.fGauss;
            obj.fValues  = cParams.fValues;
        end

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = obj.mesh.interpolation.order;
            obj.field = Field(s);
        end

        function createDiscontinousField(obj)
            s.connec = obj.meshCont.connec;
            s.type   = obj.meshCont.type;
            s.fNodes = obj.fValues';
            f = FeFunction(s);            
            fD = f.computeDiscontinousField();
            obj.fieldDisc = fD;
        end

        function createSymmetricMapCond(obj)   
            s.meshCont  = obj.meshCont;
            s.meshDisc  = obj.mesh;
            s.fieldDisc = obj.fieldDisc;
            s = SymmetricContMapCondition(s);            
            sC = s.computeCondition();
            obj.symCond = sC;            
        end        
%         
        function computeLHS(obj)
            K = obj.computeStiffnessMatrix();
            In = obj.symCond;
            K = In'*K*In;
           % M = obj.computeMassMatrix();
            I = ones(size(K,1),1);
            %eta = 0.01;
            %obj.LHS = K + eta*M;
            obj.LHS = [K,I;I',0];
        end
        
        function K = computeStiffnessMatrix(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.type         = 'StiffnessMatrix';
            s.field        = obj.field;

            lhs = LHSintegrator.create(s);
            K = lhs.compute();
        end
        
%         function M = computeMassMatrix(obj)
%             s.mesh         = obj.mesh;
%             s.globalConnec = obj.mesh.connec;
%             s.type         = 'MassMatrix';
%             s.dim          = obj.dim;
%             s.quadType     = 'QUADRATIC';
%             lhs = LHSintegrator.create(s);
%             M = lhs.compute();
%         end
        
        function computeRHS(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            s.fType     = 'Gauss';
            s.fGauss    = obj.fGauss;
            s.xGauss    = q.posgp;
            s.mesh      = obj.mesh;
            s.type      = obj.mesh.type;
            s.quadOrder = q.order;
            s.npnod     = obj.field.dim.ndofs;
            s.type      = 'ShapeDerivative';
            s.globalConnec = obj.mesh.connec;
            rhs  = RHSintegrator.create(s);
            rhsV = rhs.compute();
            In = obj.symCond;            
            rhsV = In'*rhsV;
            obj.RHS = [rhsV;0];
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
            u = u(1:end-1);
        end
        
    end
    
end