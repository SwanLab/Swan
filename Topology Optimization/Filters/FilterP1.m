classdef FilterP1 < handle

    properties (Access = private)
        mesh
     %   Poper
        M
        I
        filteredField        
        testFunction
    end

    methods (Access = public)

        function obj = FilterP1(cParams)
            obj.init(cParams);
          %  obj.createPoperator();
            obj.createMassMatrix();
            obj.createSupportMatrix();
        end

        function xReg = compute(obj,fun,quadType)
            RHS = computeRHS(obj,fun,quadType);            % computeNew
            Iki        = obj.I;
            sM         = sum(obj.M,2);
            LHS        = Iki*sM;
            P          = Iki./LHS;
            xR         = P*RHS;
            obj.filteredField.fValues(:,1) = xR;
        %    obj.filteredField.fValues(1,1,:) = xR;
            xReg = obj.filteredField;
        end
      
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.filteredField = P1Function.create(obj.mesh,1); % trial will come from outside
            obj.testFunction  = P0Function.create(obj.mesh,1);
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
% 
%         function xReg = getP1Function(obj,fun,quadType)
%             test       = P0Function.create(obj.mesh, 1);
%             int        = obj.computeRHSintegrator(quadType);
%             P          = obj.Poper.value;
%             A          = P';
%             b          = int.compute(fun,test);
%             p.fValues  = A*b;
%             p.mesh     = obj.mesh;
%             xReg       = P1Function(p);
%         end
% 
%         function xReg = getFGaussFunction(obj,fun,quadType)
%             test       = P1Function.create(obj.mesh, 1);
%             int        = obj.computeRHSintegrator(quadType);
%             P          = obj.Poper.value;
%             A          = P;
%             b          = int.compute(fun,test);
%             xR         = A*b;
%             p.fValues  = xR;
%             p.mesh     = obj.mesh;
%             xReg       = P0Function(p);
%         end

    end

end