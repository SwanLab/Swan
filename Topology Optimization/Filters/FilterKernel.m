classdef FilterKernel < handle

    properties (Access = private)
        mesh
        filteredField
        testField
    end

    properties (Access = private)
        massMatrix
        supportMatrix
        RHS
    end

    methods (Access = public)
        function obj = FilterKernel(cParams)
            obj.init(cParams);
            obj.createMassMatrix();
            obj.createSupportMatrix();
        end

        function xReg = compute(obj,fun,quadType)
            obj.computeRHS(fun,quadType);
            obj.solveFilter();
            xReg = obj.filteredField;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.filteredField = cParams.trial;
            obj.testField     = cParams.test;
        end

        function createMassMatrix(obj)
            s.type            = 'MassMatrix';
            s.mesh            = obj.mesh;
            s.test            = obj.testField;
            s.trial           = obj.testField;
            s.quadratureOrder = 'QUADRATICMASS';
            LHS               = LHSintegrator.create(s);
            obj.massMatrix    = LHS.compute();
        end      

        function createSupportMatrix(obj)
            connecTrial = obj.filteredField.computeDofConnectivity();
            connecTest  = obj.testField.computeDofConnectivity();
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
                    T      = boolean(T + incT);
                end
            end
            obj.supportMatrix = T;
        end

        function computeRHS(obj,fun,quadType)
            switch class(fun)
                case {'UnfittedFunction','UnfittedBoundaryFunction'}
                    s.mesh = fun.unfittedMesh;
                otherwise
                    s.mesh = obj.mesh;
            end
            s.type     = 'ShapeFunction';
            s.quadType = quadType;
            rhsI       = RHSintegrator.create(s);
            test       = obj.testField;
            obj.RHS    = rhsI.compute(fun,test);     
        end

        function solveFilter(obj)
            RHSi = obj.RHS;
            Iki  = obj.supportMatrix;
            Mi   = sum(obj.massMatrix);            
            Pki  = Iki./(Iki*Mi');
            xRk  = Pki*RHSi;
            obj.filteredField.fValues = xRk;
        end

    end

end