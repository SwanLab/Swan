classdef Jump < FeFunction

    properties (Access = private)
        jumpDim
        L
        R
    end

    properties (Access = private)
        uFun
        cohesiveMesh
    end

    properties (Access = public)
        fun
        nDofsElem
    end

    methods (Access = public)

        function obj = Jump(cParams)     
            obj.init(cParams);
            obj.createJumpFunction();
            obj.computeGlobalSeparationMatrix();
            obj.updateJumpValues(obj.uFun); 
        end

        function updateJumpValues(obj,uIn)          
            obj.updateRotationMatrix(uIn); % 2x2xnElems
            connJump  = obj.mesh.connec;
            nNodeJump = obj.mesh.nnodes;
            fValuesJ   = zeros(nNodeJump,obj.jumpDim); % nNodeJump x 2
            nnodesElemU = uIn.nDofsElem/uIn.ndimf;
            for n = 1:nnodesElemU 
                uNode = obj.computeDispNodes(uIn,n);
                jumpi = pagemtimes(obj.L(:,:,n*uIn.ndimf), pagemtimes(obj.R, uNode));
                jumpi = squeeze(jumpi).';
                nJump = min(n,nnodesElemU-n+1);
                fValuesJ(connJump(:,nJump),:) = fValuesJ(connJump(:,nJump),:) + jumpi;
            end
            div = [1, 2*ones(1,obj.fun.mesh.nelem-1), 1];
            fValuesJ = fValuesJ./div.';
            obj.fun.setFValues(fValuesJ);
        end

        function Bc = computeShapeFunctions(obj,xV)
            N  =  obj.fun.computeShapeFunctions(xV);  % N1(-1) N1(1); N2(-1), N2(1)
            ngauss = size(xV,2);
            nelem = obj.fun.mesh.nelem;
            Bc = zeros(obj.fun.ndimf, obj.fun.nDofsElem, ngauss, nelem);
            nnodesElemU = obj.uFun.nDofsElem/obj.uFun.ndimf;
            for i=1:obj.uFun.nDofsElem
                % N = nnode x ngauss; L = ndimf x ndimf x nelem; R = ndimf x ndimf x nelem 
                dimf = mod(i-1, obj.fun.ndimf)+1;
                nodeU    = ceil(i/obj.uFun.ndimf);
                nodeJump = min(nodeU,nnodesElemU-nodeU+1);
                Ri  = obj.R(:,dimf,:);
                Li  = obj.L(:,:,i);
                Ni  = N(nodeJump,:);
                for j = 1:ngauss % Bci = ndimf x ngauss x nelem
                    Bci(:,j,:) = Ni(j) *  pagemtimes(Li,Ri);
                end
                Bc(:,i,:,:) = Bci; %Bc = ndimf(2) x nDofElem(8) x nGauss(2) x nElem
            end
        end

        function jc = getDofConnec(obj)
            uc = obj.uFun.getDofConnec(); 
            jc = uc(obj.cohesiveMesh.listCohesiveElems,:);
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.cohesiveMesh = cParams.cohesiveMesh;
            obj.ndimf   = cParams.ndimf;
            obj.uFun    =  cParams.uFun;
            obj.jumpDim = 2;
            obj.nDofsElem = obj.uFun.nDofsElem;
            obj.mesh = obj.cohesiveMesh.mesh;
            obj.fValues = zeros(size(obj.uFun.fValues));
        end 

        function createJumpFunction(obj)
            obj.fun = LagrangianFunction.create(obj.mesh,obj.ndimf,'P1');
        end

        function computeGlobalSeparationMatrix(obj)
            obj.L = cat(3, -repmat(eye(obj.jumpDim),1,1,4), repmat(eye(obj.jumpDim),1,1,4));
        end
        
        function updateRotationMatrix(obj,uIn) 
            nCohElem     = length(obj.cohesiveMesh.listCohesiveElems);
            Rall  = zeros(obj.uFun.ndimf,obj.uFun.ndimf,nCohElem);
            for j=1:nCohElem
                e = obj.cohesiveMesh.listCohesiveElems(j);
                connecMesh = obj.cohesiveMesh.fullMesh.connec(e,:);
                coordsMesh = obj.cohesiveMesh.fullMesh.coord(connecMesh',:);
                disp  = uIn.fValues(connecMesh',:);
                    deformedCoords = coordsMesh + disp;
                Re = obj.computeElementalRotationMatrix(deformedCoords);
                Rall(:,:,j) = Re;
                obj.R = Rall;
            end
        end

        function uVals = computeDispNodes(obj,uIn,n)
            conn  = obj.cohesiveMesh.fullMesh.connec(obj.cohesiveMesh.listCohesiveElems,:);
            nodes = conn(:,n);
            uVals = uIn.fValues(nodes,:);
            uVals = reshape(uVals.',[2 1 size(uVals,1)]);
        end

        function Re = computeElementalRotationMatrix(obj,deformedCoords)
                midPoints= 0.5*[deformedCoords(1,1)+deformedCoords(4,1),deformedCoords(1,2)+deformedCoords(4,2);
                             deformedCoords(2,1)+deformedCoords(3,1),deformedCoords(2,2)+deformedCoords(3,2)];
                m = [midPoints(2,1)-midPoints(1,1),midPoints(2,2)-midPoints(1,2)];
                mx = m(1); my = m(2);   
                Re = [mx, -my; my, mx] / sqrt(mx^2 + my^2);
        end
    end

    methods(Access = protected)

        function fV = evaluateNew(obj,xV)
        fV = obj.fun.evaluate(xV);
        end       

    end
end
