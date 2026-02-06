classdef CohesiveSeparationComputer < handle
    properties (Access = public)
        cohesiveMesh
        globalSeparations % nCohElem x nMidNodesPerElem(2) x nSeparations(2)
        lagrangianSeparation
        
        fun

        ndimf
        order

        L
    end

    methods (Access = public)

        function obj = CohesiveSeparationComputer(cParams)

            obj.init(cParams);
            obj.createJumpFunction %(initialize the fun with zeros
            obj.L = obj.computeGlobalSeparationMatrix();
        end

        function compute(obj,uIn)
            uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);

            R = obj.rotationMatrix(uIn);

            funOutVec = obj.L * R * uInVec; % ndofDisp x 1

            funOut = funOutVec;

        end

        function fV = evaluate(obj,xV)
            fV = obj.fun.evaluate(xV);
        end

        % validate by setting down face displacement to: 1) x=-1, 2) y=1, 3) x=[-1:-4], 4) y[-1:-4]  
        
        function globalSeps = ComputeGlobalSeparations(obj)
            %%%% uInVec = reshape(uIn.fValues',[uIn.nDofs 1]); 
            %%%% uOut = reshape(uOutVec,[flip(size(uIn.fValues))])';
            connec = obj.cohesiveMesh.mesh.connec;   % (nElem × 4)
            coords = obj.cohesiveMesh.mesh.coord;    % (nNodes × ndim)
        
            L = [-1 0 0 1;
                  0 -1 1 0];                         % (2×4)
        
            nelem = length(obj.cohesiveMesh.listCohesiveElems);
            globalSeps = zeros(nelem,2,2);

            for i = 1:nelem
                e = obj.cohesiveMesh.listCohesiveElems(i);
                R = obj.rotationMatrix(i);


                
                Ue = obj.u.fValues(connec(e,:)',:);
                globalSeps(i,:,:) = L*Ue*R';
            end
            
        end

   end


    methods (Access = private)
        function init(obj,cParams)
            obj.cohesiveMesh = cParams.cohesiveMesh;
            obj.ndimf = cParams.ndimf;
            % obj.order = cParams.order;
        end

        function R = rotationMatrix(obj,uIn)

            nDofsU = uIn.nDofs;
            R = sparse(nDofsU,nDofsU);
            nCohElem = length(obj.cohesiveMesh.listCohesiveElems);

            for j=1:nCohElem
                e = obj.cohesiveMesh.listCohesiveElems(j);
                connecSub = obj.cohesiveMesh.subMesh.connec(j,:);

                coords = obj.cohesiveMesh.subMesh.coord(connecSub',:)+uIn.fValues(connecSub',:);

                ndim = obj.cohesiveMesh.mesh.ndim;
                
                connecFull = obj.cohesiveMesh.mesh.connec(e,:);
                dofsU =  reshape((connecFull(:)-1)*ndim + (1:ndim), 1, []);

                m = [coords(2,1)-coords(1,1),coords(2,2)-coords(1,2)];
                    mx = m(1);
                    my = m(2);
                    
                Re = [mx, -my; my, mx] / sqrt(mx^2 + my^2);
                ReBig = kron(eye(4), Re); %generalitzar el 4

                R(dofsU,dofsU) = ReBig;

            end

        end

        function lagrangianSeparation = ComputeLagrangianSeparation(obj)

            fValues = globalSeps';
            obj.fun.setFValues(fValues);

        end

        function createJumpFunction(obj)
            obj.fun = LagrangianFunction.create(obj.cohesiveMesh.subMesh,2,'P1D');
        end


        function L = computeGlobalSeparationMatrix(obj)
            % L -- ndofJump x ndofu
            nCohElem = length(obj.cohesiveMesh.listCohesiveElems);
        
            nJumpPerElem = 4;

            nDofU = obj.cohesiveMesh.mesh.nnodes * obj.cohesiveMesh.mesh.ndim;
            nDofJump = 4* nCohElem; % arreglar
      
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
                dofsJump = 2*(j-1)*ones(nJumpPerElem,1).'+(1:nJumpPerElem);

                L(dofsJump,dofsU) = Le;

            end







        end


    end


end
