classdef FilterP1 < handle

    properties (Access = private)
        mesh
       Poper
        M
        M2
        I
        filteredField        
        testFunction
    end

    methods (Access = public)

        function obj = FilterP1(cParams)
            obj.init(cParams);
           obj.createPoperator();
            obj.createMassMatrix();
            obj.createMassMatrix2();
            obj.createSupportMatrix();
        end

        function xReg = compute(obj,fun,quadType)
            RHS  = obj.computeRHS(fun,quadType);
            Iki  = obj.I;
            %             It  = ones(size(obj.M,2),1);
            %             LHS = Iki*obj.M*It;
            %             xR  = diag(LHS)\(Iki*RHS);
            It2  = ones(size(obj.M2,2),1);
            LHS2 = Iki*obj.M2*It2;
            xR   = diag(LHS2)\(Iki*RHS);
            obj.filteredField.fValues(:) = xR;
            xReg = obj.filteredField;
        end

        function xReg = compute2(obj,fun,quadType) % computeWorking DELETE ASAP
            switch class(fun)
                case 'FGaussDiscontinuousFunction'
                    xReg = obj.getP1Function(fun,quadType);
                otherwise
                    xReg = obj.getFGaussFunction(fun,quadType);
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.filteredField = cParams.trial;
            obj.testFunction  = cParams.test;
        end

        function createPoperator(obj)
            s.mesh    = obj.mesh;
            obj.Poper = Poperator(s);
        end

        function createMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.testFunction;
            s.trial = obj.testFunction;
            s.quadratureOrder = 'QUADRATICMASS';
            LHS   = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end

        function createMassMatrix2(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.testFunction;
            s.trial = obj.filteredField;
            s.quadratureOrder = 'QUADRATICMASS';
            LHS   = LHSintegrator.create(s);
            obj.M2 = LHS.compute();
        end        

        function createSupportMatrix(obj)
            connecTrial = obj.filteredField.computeDofConnectivity();
            connecTest  = obj.testFunction.computeDofConnectivity();
            nDofsP1     = max(connecTest, [], 'all');
            nDofsField  = max(connecTrial, [], 'all');
            nDofElemP1  = size(connecTest,1);
            nDofElemF   = size(connecTrial,1);
            T = sparse(nDofsField,nDofsP1);
            for kDof = 1:nDofElemF
                for iDof = 1:nDofElemP1
                    dofsF  = connecTrial(kDof,:);
                    dofsP1 = connecTest(iDof,:);
                    Iv     = ones(obj.mesh.nelem,1);
                    incT   = sparse(dofsF,dofsP1,Iv,nDofsField,nDofsP1);
                    T      = T + incT;
                end
            end
            obj.I = T;
        end

        function RHS = computeRHS(obj,fun,quadType)
            s.type     = 'ShapeFunction';
            s.quadType = quadType;
            s.mesh     = obj.mesh;
            rhsI       = RHSintegrator.create(s);            
            test       = obj.testFunction;
            RHS  = rhsI.compute(fun,test);     
        end

        function xReg = getP1Function(obj,fun,quadType)
            b        = obj.computeRHS(fun,quadType);
            P          = obj.Poper.value;
            A          = P';
            p.fValues  = A*b;
            p.mesh     = obj.mesh;
            xReg       = P1Function(p);
        end

        function xReg = getFGaussFunction(obj,fun,quadType)
            b        = obj.computeRHS(fun,quadType);
            P          = obj.Poper.value;
            A          = P;
            xR         = A*b;
            p.fValues  = xR;
            p.mesh     = obj.mesh;
            xReg       = P0Function(p);
        end

    end

end