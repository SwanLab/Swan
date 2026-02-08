classdef CohesiveSeparationComputer < handle
    properties (Access = public)
        cohesiveMesh
        globalSeparations % nCohElem x nMidNodesPerElem(2) x nSeparations(2)
        lagrangianSeparation
        
        fun
        effectiveFun

        jumpDim

        ndimf
        order

        L
    end

    methods (Access = public)

        function obj = CohesiveSeparationComputer(cParams)
            cParams.jumpDim = 2;
            obj.init(cParams);
            obj.createJumpFunction();
            obj.L = obj.computeGlobalSeparationMatrix();
        end

        function compute(obj,uIn)
            uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);

            R = obj.rotationMatrix(uIn);

            funOutVec = obj.L * R * uInVec; % ndofDisp x 1

            funOut = reshape(funOutVec,[obj.jumpDim obj.fun.nDofs/obj.jumpDim]);

            obj.updateJumps(funOut);


        end

        function fV = evaluate(obj,xV)
            fV = obj.fun.evaluate(xV);
        end


        function computeEffectiveSeps(obj,uIn)
            obj.compute(uIn);
            fValues = obj.fun.fValues;
            effectiveFValues = vecnorm(fValues,2,2);
            obj.effectiveFun = LagrangianFunction.create(obj.cohesiveMesh.subMesh,obj.jumpDim,'P1D');
            obj.effectiveFun.setFValues(effectiveFValues);
        end

   end


    methods (Access = private)
        function init(obj,cParams)
            obj.cohesiveMesh = cParams.cohesiveMesh;
            obj.ndimf = cParams.ndimf;
            obj.jumpDim = cParams.jumpDim;
            % obj.order = cParams.order;
            
        end

        function R = rotationMatrix(obj,uIn)

            nDofsU = uIn.nDofs;
            R = sparse(nDofsU,nDofsU);
            nCohElem = length(obj.cohesiveMesh.listCohesiveElems);
            nJumpPerElem = obj.jumpDim * obj.cohesiveMesh.subMesh.nnodeElem;
            
            for j=1:nCohElem
                e = obj.cohesiveMesh.listCohesiveElems(j);
                
                connecMesh = obj.cohesiveMesh.mesh.connec(e,:);
                coordsMesh = obj.cohesiveMesh.mesh.coord(connecMesh',:);

                dofsU = reshape(((connecMesh(:)-1)*obj.ndimf + (1:obj.ndimf)).', 1, []);
               
                disp = uIn.fValues(connecMesh',:);

                deformedCoords = coordsMesh + disp;

                midPoints = 0.5*[deformedCoords(1,1)+deformedCoords(4,1),deformedCoords(1,2)+deformedCoords(4,2);
                             deformedCoords(2,1)+deformedCoords(3,1),deformedCoords(2,2)+deformedCoords(3,2)];

                m = [midPoints(2,1)-midPoints(1,1),midPoints(2,2)-midPoints(1,2)];
                    mx = m(1);
                    my = m(2);
                    
                Re = [mx, -my; my, mx] / sqrt(mx^2 + my^2);
                ReBig = kron(eye(nJumpPerElem), Re);

                R(dofsU,dofsU) = ReBig;
            end

        end

        function updateJumps(obj,funOut)

            fValues = funOut';
            obj.fun.setFValues(fValues);

        end


        function createJumpFunction(obj)
            obj.fun = LagrangianFunction.create(obj.cohesiveMesh.subMesh,obj.jumpDim,'P1D');
        end


        function L = computeGlobalSeparationMatrix(obj)
            % L -- ndofJump x ndofu
            nCohElem = length(obj.cohesiveMesh.listCohesiveElems);
        
            nJumpPerElem = obj.jumpDim * obj.cohesiveMesh.subMesh.nnodeElem;

            nDofU = obj.cohesiveMesh.mesh.nnodes * obj.cohesiveMesh.mesh.ndim;
            nDofJump = nJumpPerElem* nCohElem;
      
            ndim = obj.cohesiveMesh.mesh.ndim;

            L = sparse(nDofJump, nDofU);

            Le = [-1,0,0,0,0,0,1,0;
                   0,-1,0,0,0,0,0,1;
                   0,0,-1,0,1,0,0,0;
                   0,0,0,-1,0,1,0,0];
     
            for j = 1:nCohElem
                elem = obj.cohesiveMesh.listCohesiveElems(j);
                connec = obj.cohesiveMesh.mesh.connec(elem,:);

                dofsU = reshape(((connec(:)-1)*ndim + (1:ndim)).', 1, []);
                dofsJump = 4*(j-1)*ones(nJumpPerElem,1).'+(1:nJumpPerElem);

                L(dofsJump,dofsU) = Le;
            end

        end
        
    end


end
