classdef ClementInterpolator < handle

    properties (Access = private)
        mesh
        trial
        test
    end

    properties (Access = private)
        massMatrix
        supportMatrix
        localMeshes
        local2GlobalDofs
        RHS
        IElems
    end

    methods (Access = public)
        function obj = ClementInterpolator(cParams)
            obj.init(cParams);
            obj.createNeighborElementsMatrix();
            obj.createSupportMatrix();
            obj.createLocalMeshes();
            obj.createLocalMassMatrices();
        end

        function xReg = compute(obj,fun,quadType)
            obj.computeLocalRHS(fun,quadType);
            obj.solve();
            xReg = obj.trial;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.trial = P1Function.create(obj.mesh,1);
            obj.test  = P1Function.create(obj.mesh,1);
        end

        function createSupportMatrix(obj)
            connecTrial = obj.trial.computeDofConnectivity();
            p0          = P0Function.create(obj.mesh,1);
            connecElem  = p0.computeDofConnectivity;
            nDofsP0     = obj.mesh.nelem;
            nDofsField  = max(connecTrial, [], 'all');
            nDofElemP0  = size(connecElem,1);
            nDofTrial   = size(connecTrial,1);
            T           = sparse(nDofsField,nDofsP0);
            for kDof = 1:nDofTrial
                for iDof = 1:nDofElemP0
                    dofsF  = connecTrial(kDof,:);
                    dofsP0 = connecElem(iDof,:);
                    Iv     = ones(nDofsP0,1);
                    incT   = sparse(dofsF,dofsP0,Iv,nDofsField,nDofsP0);
                    T      = boolean(T + incT);
                end
            end
            obj.supportMatrix = min(T,1); % obj.IElems !!! generalize
        end

        function createLocalMeshes(obj)
            connecElem  = obj.mesh.connec;
            nelem       = obj.mesh.nelem;
            nnode       = obj.mesh.nnodes;
            nnodeElem   = size(connecElem,2);
            T           = obj.supportMatrix;
            Tsort       = reshape(full(T'),[nelem,1,nnode]);
            localConnec = connecElem.*Tsort;

            % Pending to vectorize-----------

            for i = 1:nnode
                lC       = localConnec(:,:,i);
                lC       = nonzeros(lC');
                locnElem = length(lC)/nnodeElem;
                lC       = reshape(lC,[nnodeElem,locnElem])';
                s.connec = lC;
                s.coord  = obj.mesh.coord;
                msh      = Mesh(s);
                msh      = msh.computeCanonicalMesh();
                locCoord = msh.coord;
                globNode = obj.mesh.coord(i,:);
                obj.local2GlobalDofs(i) = find(ismember(locCoord,globNode,'rows'));
                obj.localMeshes{i}      = msh;
            end
            % --------------------
        end

        function createLocalMassMatrices(obj)
            s.type            = 'MassMatrix';
            s.quadratureOrder = 'QUADRATICMASS';
            for i = 1:obj.mesh.nnodes
                msh               = obj.localMeshes{i};
                s.test            = P1Function.create(msh,1);
                s.trial           = P1Function.create(msh,1);
                s.mesh            = msh;
                LHS               = LHSintegrator.create(s);
                obj.massMatrix{i} = LHS.compute();
            end
        end

        function computeLocalRHS(obj,fun,quadType)
            s.type     = 'ShapeFunction';
            s.quadType = quadType;
            for i = 1:obj.mesh.nnodes
                msh = obj.localMeshes{i};
                fun.updateMesh(msh);
                switch class(fun)
                    case {'UnfittedFunction','UnfittedBoundaryFunction'}
                        s.mesh = fun.unfittedMesh;
                    otherwise
                        s.mesh = msh;
                end
                rhsI       = RHSintegrator.create(s);
                obj.RHS{i} = rhsI.compute(fun,P1Function.create(msh,1));
            end
        end

        function solve(obj)
            locM   = obj.massMatrix;
            locRHS = obj.RHS;
            fVal   = zeros(size(obj.trial.fValues));
            for i = 1:obj.mesh.nnodes
                Mi      = locM{i};
                Mi      = obj.lumpMatrix(Mi);
                RHSi    = locRHS{i};
                xi      = Mi\RHSi;
                fVal(i) = xi(obj.local2GlobalDofs(i));
            end

            obj.trial.fValues = fVal;
        end

        function Al = lumpMatrix(obj,A)
            I  = ones(size(A,2),1);
            Al = diag(A*I);
        end

        function createNeighborElementsMatrix(obj)
            connec       = obj.mesh.connec;
            nodes(1,1,:) = 1:obj.mesh.nnodes;
            neig         = sum(connec==nodes,2);
            neig         = squeeze(neig);
            T            = neig*neig';
            obj.IElems   = min(T,1);
        end

    end

end