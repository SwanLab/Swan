classdef FilterKernelGeneral < handle

    properties (Access = private)
        mesh
        filteredField
        testField
    end

    properties (Access = private)
        massMatrix
        IElems
        supportMatrix
        RHS
    end

    methods (Access = public)
        function obj = FilterKernelGeneral(cParams)
            obj.init(cParams);
            obj.createMassMatrix();
            obj.createNeighborElementsMatrix();
            obj.createGlobalSupportMatrix();
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
            s.trial           = obj.filteredField;
            s.quadratureOrder = 'QUADRATICMASS';
            LHS               = LHSintegrator.create(s);
            obj.massMatrix    = LHS.compute();
        end 

        function createNeighborElementsMatrix(obj)
            connec       = obj.mesh.connec;
            nodes(1,1,:) = 1:obj.mesh.nnodes;
            neig         = sum(connec==nodes,2); % Be careful! Maybe shared edges, not shared nodes
            neig         = squeeze(neig);
            T            = neig*neig';
            obj.IElems   = min(T,1);
        end

        function createGlobalSupportMatrix(obj)
            elems             = 1:obj.mesh.nelem;
            locSM             = obj.createLocalSupportMatrix(elems);
            nElem             = obj.mesh.nelem;
            IEl               = reshape(obj.IElems,[nElem,1,nElem]);
            neigsElem         = sum(squeeze(IEl),2)-1;
            maxNeig           = max(neigsElem);
            for i = 1:maxNeig
                iNeig(1,1,:) = neigsElem>=i;
                IEl.*iNeig % here I stayed
            end
            obj.supportMatrix = min(locSM,1);
        end

        function locSM = createLocalSupportMatrix(obj,elems)
            connecTrial = obj.filteredField.computeDofConnectivity();
            connecTest  = obj.testField.computeDofConnectivity();
            nDofsP1     = max(connecTest, [], 'all');
            nDofsField  = max(connecTrial, [], 'all');
            nDofElemP1  = size(connecTest,1);
            nDofElemF   = size(connecTrial,1);
            T = sparse(nDofsField,nDofsP1);
            for kDof = 1:nDofElemF
                for iDof = 1:nDofElemP1
                    dofsF  = connecTrial(kDof,elems);
                    dofsP1 = connecTest(iDof,elems);
                    Iv     = ones(obj.mesh.nelem,1);
                    incT   = sparse(dofsF,dofsP1,Iv,nDofsField,nDofsP1);
                    T      = T + incT;
                end
            end
            locSM = min(T,1);
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
            rhs = obj.RHS;
            Iki = obj.supportMatrix;
            M   = obj.massMatrix;
            LHS = Iki*M;
            LHS = obj.lumpMatrix(LHS);
            xRk = (Iki*rhs)./LHS;
            obj.filteredField.fValues = xRk;
        end

    end

    methods (Static, Access = private)

        function Al = lumpMatrix(A)
            I  = ones(size(A,2),1);
            Al = A*I;
        end

    end

end