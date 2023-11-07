classdef FilterP1 < handle

    properties (Access = private)
        mesh
        filteredField
        Poper
        M
        I
    end

    methods (Access = public)

        function obj = FilterP1(cParams)
            obj.init(cParams);
            obj.createPoperator();
            obj.createMassMatrix();
            obj.createSupportMatrix();
        end

        function xReg = computeNew(obj,fun,quadType) % computeNew
            test       = P1Function.create(obj.mesh, 1);
            int        = obj.computeRHSintegrator(quadType);
            Iki        = obj.I;
            sM         = sum(obj.M,2);
            den        = Iki*sM;
            P          = Iki./den;
            A          = P;
            b          = int.compute(fun,test);
            xR         = A*b;
            p.fValues  = xR;
            p.mesh     = obj.mesh;
            xReg       = P1Function(p);
        end

        function xReg = compute(obj,fun,quadType) % computeWorking DELETE ASAP
            switch class(fun)
                case 'P1Function'
                    xReg = obj.getFGaussFunction(fun,quadType);
                case 'FGaussDiscontinuousFunction'
                    xReg = obj.getP1Function(fun,quadType);
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.filteredField = P1Function.create(obj.mesh,1); % trial will come from outside
        end

        function createPoperator(obj)
            s.mesh    = obj.mesh;
            obj.Poper = Poperator(s);
        end

        function createMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.quadratureOrder = 'QUADRATICMASS';
            LHS   = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end

%         function createNeighborElementsMatrix(obj)
%             nelem     = obj.mesh.nelem;
%             nodesElem = obj.mesh.nnodeElem;
%             connec    = obj.mesh.connec;
%             T = zeros(nelem,nelem);
%             for ielem = 1:nelem
%                 for inode=1:nodesElem
%                     node = connec(ielem,inode);
%                     neigElems = connec==node;
%                     neigElems = any(neigElems');
%                     T(ielem,neigElems) = 1;
%                 end
%             end
%             obj.IElems = T;
%         end

        function createSupportMatrix(obj)
            connecField = obj.filteredField.computeDofConnectivity;
            connecP1    = obj.mesh.connec';
            nelem       = obj.mesh.nelem;
            nnodes      = obj.mesh.nnodes;
            nDofs       = max(connecField, [], 'all');
            nodesElem   = obj.mesh.nnodeElem;
            nDofElem    = size(connecField,1);
            T           = zeros(nDofs,nnodes);
            for ielem = 1:nelem
                for kdof = 1:nDofElem
                    for inode = 1:nodesElem
                        dofs  = connecField(kdof,ielem);
                        nodes = connecP1(inode,ielem);
                        T(dofs,nodes) = 1;
                    end
                end
            end
            obj.I = T;
        end

        function rhs = computeRHSintegrator(obj,quadType)
            s.type     = 'ShapeFunction';
            s.quadType = quadType;
            s.mesh     = obj.mesh;
            rhs        = RHSintegrator.create(s);
        end

        function xReg = getP1Function(obj,fun,quadType)
            test       = P0Function.create(obj.mesh, 1);
            int        = obj.computeRHSintegrator(quadType);
            P          = obj.Poper.value;
            A          = P';
            b          = int.compute(fun,test);
            p.fValues  = A*b;
            p.mesh     = obj.mesh;
            xReg       = P1Function(p);
        end

        function xReg = getFGaussFunction(obj,fun,quadType)
            test       = P1Function.create(obj.mesh, 1);
            int        = obj.computeRHSintegrator(quadType);
            P          = obj.Poper.value;
            A          = P;
            b          = int.compute(fun,test);
            xR         = A*b;
            p.fValues  = xR;
            p.mesh     = obj.mesh;
            xReg       = P0Function(p);
        end

    end

end