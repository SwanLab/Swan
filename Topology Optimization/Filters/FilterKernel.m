classdef FilterKernel < handle

    properties (Access = private)
        mesh
        filteredField
        testField
        approach
    end

    properties (Access = private)
        massMatrix
        M2
        M3
        M4
        supportMatrix
        RHS
        RHS2
    end

    methods (Access = public)
        function obj = FilterKernel(cParams)
            obj.init(cParams);
            obj.createMassMatrix();
            obj.createSupportMatrix();
        end

        function xReg = compute(obj,fun,quadType)
            obj.computeRHS(fun,quadType);
            obj.computeRHS2(fun,quadType);
            obj.solveFilter();
            xReg = obj.filteredField;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.filteredField = cParams.trial;
            obj.testField     = cParams.test;
            if not(isfield(cParams,'approach'))
                obj.approach      = 'B';
            end
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

        function createMassMatrix2(obj)
            s.type            = 'MassMatrix';
            s.mesh            = obj.mesh;
            s.test            = obj.testField;
            s.trial           = obj.filteredField;
            s.quadratureOrder = 'QUADRATICMASS';
            LHS               = LHSintegrator.create(s);
            obj.M2    = LHS.compute();
        end 

        function createMassMatrix3(obj)
            s.type            = 'MassMatrix';
            s.mesh            = obj.mesh;
            s.test            = obj.filteredField;
            s.trial           = obj.filteredField;
            s.quadratureOrder = 'QUADRATICMASS';
            LHS               = LHSintegrator.create(s);
            obj.M3    = LHS.compute();
        end   

        function createMassMatrix4(obj)
            s.type            = 'MassMatrix';
            s.mesh            = obj.mesh;
            s.test            = obj.filteredField;
            s.trial           = obj.testField;
            s.quadratureOrder = 'QUADRATICMASS';
            LHS               = LHSintegrator.create(s);
            obj.M4    = LHS.compute();
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

        function computeRHS2(obj,fun,quadType)
            switch class(fun)
                case {'UnfittedFunction','UnfittedBoundaryFunction'}
                    s.mesh = obj.mesh;
                otherwise
                    s.mesh = obj.mesh;
            end
            s.type     = 'ShapeFunction';
            s.quadType = quadType;
            rhsI       = RHSintegrator.create(s);
            test       = obj.filteredField;
            obj.RHS2    = rhsI.compute(fun,test);
        end

        function solveFilter(obj)
            obj.createMassMatrix2();
            obj.createMassMatrix3();
            obj.createMassMatrix4();

            RHSi = obj.RHS;
            rhs2 = obj.RHS2;

            Iki  = obj.supportMatrix;

            switch (obj.approach)

                case {'B'}

                    M2   = obj.M2;
                    %I1     = ones(size(M,2),1);
                    %I2     = ones(size(M2,2),1);
                    %LHSp  = Iki*M;
                    LHS = Iki*M2;
                    %Mi   = sum(M,2);
                    %Mi  = M*I1;
                    %LHS  = diag(LHSp*I1);
                    %LHS2 = diag(LHSp2*I2);
                    LHS  = obj.lumpMatrix(LHS);
                    %norm(LHS(:)-LHS2(:))/norm(LHS(:))
                    xRk  = (Iki*RHSi)./LHS;

                case {'A'}

                    M3    = obj.M3;
                    LHS2 = Iki'*M3*Iki;
                    LHS2 = obj.lumpMatrix(LHS2);
                    xRk  = (Iki'*rhs2)./LHS2;
                    xRk = (Iki)*xRk;

                    %norm(xRk - xRk2)/norm(xRk)
                    %  xRk2  = Iki'*xRk2;

                case {'C'}

                    M2   = obj.M2;
                    %I1     = ones(size(M,2),1);
                    %I2     = ones(size(M2,2),1);
                    %LHSp  = Iki*M;
                    LHS = Iki*M2;
                    %Mi   = sum(M,2);
                    %Mi  = M*I1;
                    %LHS  = diag(LHSp*I1);
                    %LHS2 = diag(LHSp2*I2);
                    LHS  = obj.lumpMatrix(LHS);
                    %norm(LHS(:)-LHS2(:))/norm(LHS(:))
                    P    = (Iki)./LHS;
                    xRk  = P'*rhs2;
                    obj.filteredField = P1Function.create(obj.mesh,1);

            end



            obj.filteredField.fValues = xRk;
        end

        function Al = lumpMatrix(obj,A)
            I  = ones(size(A,2),1);
            Al = A*I;
        end

    end

end