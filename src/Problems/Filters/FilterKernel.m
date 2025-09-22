classdef FilterKernel < handle

    properties (Access = private)
        mesh
        trial
        test
        nLevels
    end

    properties (Access = private)
        massMatrix
        IElems
        supportMatrix
        RHS
    end

    methods (Access = public)
        function obj = FilterKernel(cParams)
            obj.init(cParams);
            obj.createMassMatrix();
            %obj.createNeighborElementsMatrix();
            obj.createGlobalSupportMatrix();
        end

        function xReg = compute(obj,fun,quadType)
            xReg = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);
            obj.computeRHS(fun,quadType);
            obj.solveFilter();
            xReg.setFValues(full(obj.trial.fValues));
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            cParams.feFunType = class(cParams.trial);
            cParams.ndimf     = 1;
            obj.trial         = LagrangianFunction.create(cParams.mesh, 1, cParams.trial.order);
            obj.test          = cParams.test;
            obj.mesh          = cParams.mesh;
            obj.nLevels       = 1; % Must be defined before entering the class
        end 

        function createMassMatrix(obj)
            f = @(v,u) DP(v,u);               
            obj.massMatrix = IntegrateLHS(f,obj.test,obj.trial,obj.mesh,2);
        end 

        function createNeighborElementsMatrix(obj)
            connec       = obj.mesh.connec;
            nodes(1,1,:) = 1:obj.mesh.nnodes;
            neig         = sum(connec==nodes,2);
            neig         = squeeze(neig);
            T            = neig*neig';
            obj.IElems   = min(T,1);
        end

        function createGlobalSupportMatrix(obj)
            elems = 1:obj.mesh.nelem;
            locSM = obj.createLocalSupportMatrix(elems);
            if obj.nLevels == 1 % Provisional solution before solving issue
                obj.supportMatrix = locSM;
            else % not generalized:
                IEl               = obj.IElems;
                p                 = obj.nLevels;
                patchMat          = IEl^(p-1);
                obj.supportMatrix = min(locSM*patchMat,1);
            end
        end

        function locSM = createLocalSupportMatrix(obj,elems)
            connecTrial  = obj.trial.getDofConnec()';
            connecTest   = obj.test.getDofConnec()';
            nDofsTest    = max(max(connecTest));%obj.test.nDofs;
            nDofsField   = max(max(connecTrial));%obj.trial.nDofs;
            nDofElemTest = size(connecTest,1);
            nDofElemF    = size(connecTrial,1);
            T            = sparse(nDofsField,nDofsTest);
            for kDof = 1:nDofElemF
                for iDof = 1:nDofElemTest
                    dofsF    = connecTrial(kDof,elems);
                    dofsTest = connecTest(iDof,elems);
                    Iv       = ones(obj.mesh.nelem,1);
                    incT     = sparse(dofsF,dofsTest,Iv,nDofsField,nDofsTest);
                    T        = T + incT;
                end
            end
            locSM = min(T,1);
        end

        function computeRHS(obj,fun,quadType)
            switch class(fun)
                case {'UnfittedFunction','UnfittedBoundaryFunction'}
                    s.mesh = fun.unfittedMesh;
                    s.type = 'Unfitted';
                    s.quadType = quadType;
                    int        = RHSIntegrator.create(s);
                    obj.RHS    = int.compute(fun,obj.test);
                otherwise
                    f = @(v) DP(v,fun);
                    obj.RHS = IntegrateRHS(f,obj.test,obj.test.mesh,quadType);   
            end            
        end

        function solveFilter(obj)
            rhs = obj.RHS;
            Iki = obj.supportMatrix;
            M   = obj.massMatrix;
            LHS = Iki*M;
            LHS = obj.lumpMatrix(LHS);
            xRk = (Iki*rhs)./LHS;
            obj.trial.setFValues(xRk);
        end

    end

    methods (Static, Access = private)

        function Al = lumpMatrix(A)
            I  = ones(size(A,2),1);
            Al = A*I;
        end

    end

end