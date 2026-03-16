classdef Jump < handle

    properties (Access = private)
        jumpDim
        L
    end

    properties (Access = private)
        uFun
        cohesiveMesh     
        ndimf 
    end

    properties (Access = public)
        jumpFun
    end

    methods (Access = public)

        function obj = Jump(cParams)     
            obj.init(cParams);
            obj.createJumpFunction();
            obj.computeGlobalSeparationMatrix();
            obj.updateJumpValues(obj.uFun); 
        end
        function updateJumpValues(obj,uIn)          
            R = obj.computeRotationMatrix(uIn); % 2x2xnElems
            connJump  = obj.cohesiveMesh.lineMesh.connec;
            nNodeJump = obj.cohesiveMesh.lineMesh.nnodes;

            % fValues   = zeros(nNodeJump,obj.jumpDim); % nNode x 2
            % 
            % for i = 1:uIn.nDofsElem
            %     uVals = obj.computeDispDofs(uIn, i); %2X1Xnelem
            %     uVals = uVals(:,:,obj.cohesiveMesh.listCohesiveElems);
            % 
            %     jumpi = pagemtimes(obj.L(:,:,i), pagemtimes(R, uVals)); 
            %     jumpi = squeeze(jumpi).';
            % 
            %     localNode  = 1 + ~(i<=2 || i>=7);
            %     globalNode = connJump(:,localNode);
            %     fValues(globalNode,:) = fValues(:,:) + jumpi;
            % end
            % obj.jumpFun.setFValues(fValues);

            fValues   = zeros(nNodeJump,obj.jumpDim); % nNode x 2
            nnodesElemU = uIn.nDofsElem/uIn.ndimf;
            for n = 1:nnodesElemU 
                uNode = computeDispNodes(obj,uIn,n);
                

                jumpi = pagemtimes(obj.L(:,:,n*uIn.ndimf), pagemtimes(R, uNode));
                jumpi = squeeze(jumpi).';
                
                nJump = min(n,nnodesElemU-n+1);

                fValues(connJump(:,nJump),:) = fValues(connJump(:,nJump),:) + jumpi;

            end

            % div = [1, 2*ones(1,nelem-2), 1];
            % fValues = fValues./div;

        end

        function fV = evaluate(obj,xV)
            fV = obj.jumpFun.evaluate(xV);
        end
        
        
        function N = computeShapeFunctions(obj,xV) % pel test % 2x8xnelem
            R  =  obj.computeRotationMatrix(obj.uFun); % 2x2xnElem
            N  =  obj.jumpFun.computeShapeFunctions(xV);  % N1(-1) N1(1); N2(-1), N2(1)
            ngauss = size(xV,2);
            nelem = obj.jumpFun.mesh.nelem;
            Bc = zeros(obj.jumpFun.ndimf, obj.jumpFun.nDofsElem, ngauss, nelem);
            nnodesElemU = obj.uFun.nDofsElem/obj.uFun.ndimf;
            for i=1:obj.uFun.nDofsElem
                % N = nnode x ngauss
                % R = ndimf x ndimf x nelem 
                % L = ndimf x ndimf x nelem
                dimf = mod(i-1, obj.jumpFun.ndimf)+1;
                nodeU    = ceil(i/obj.uFun.ndimf);
                nodeJump = min(nodeU,nnodesElemU-nodeU+1);
                Ri  = R(:,dimf,:);
                Li  = obj.L(:,:,i);
                Ni  = N(nodeJump,:);
                for j = 1:ngauss % Bci = ndimf x ngauss x nelem
                    Bci(:,j,:) = Ni(j) *  pagemtimes(Li,Ri);
                end
                Bc(:,i,:,:) = Bci;
            end
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.cohesiveMesh = cParams.cohesiveMesh;
            obj.ndimf = cParams.ndimf;
            obj.uFun  =  cParams.uFun;
            obj.jumpDim = 2;
        end

        function createJumpFunction(obj)
            obj.jumpFun = LagrangianFunction.create(obj.cohesiveMesh.lineMesh,obj.ndimf,'P1');
        end

        function computeGlobalSeparationMatrix(obj)
            obj.L =cat(3, -repmat(eye(2),1,1,4), repmat(eye(2),1,1,4));
        end
        
        function Rall = computeRotationMatrix(obj,uIn) 
            nCohElem     = length(obj.cohesiveMesh.listCohesiveElems);
            Rall  = zeros(2,2,nCohElem);

            for j=1:nCohElem
                e = obj.cohesiveMesh.listCohesiveElems(j);
                connecMesh = obj.cohesiveMesh.mesh.connec(e,:);
                coordsMesh = obj.cohesiveMesh.mesh.coord(connecMesh',:);
                disp  = uIn.fValues(connecMesh',:);
                    deformedCoords = coordsMesh + disp;
                Re = obj.elementalRotationMatrix(deformedCoords);
                Rall(:,:,j) = Re;
            end
        end
        
        function uValsfinal = computeDispDofs(obj,uIn,i) % MALAMENT!!!!
            dim = mod(i-1, uIn.ndimf) + 1;
            conn      = obj.cohesiveMesh.mesh.connec; %nelemxnnode
            node      = ceil(i/uIn.ndimf); % 1-4
            nodes     = conn(:,node); %nelemx1
           
            uVals = uIn.fValues(nodes,:);
            uVals = uVals(:,dim);
            n = length(uVals);
            v = zeros(2,1,n);
            v(sub2ind([2 n], mod((1:n)-1,2)+1, 1:n)) = uVals;
            uValsfinal = v;
        end
    

        function uVals = computeDispNodes(obj,uIn,n)
            conn = obj.cohesiveMesh.mesh.connec(obj.cohesiveMesh.listCohesiveElems,:);
            nodes = conn(:,n);
            uVals = uIn.fValues(nodes,:);
            uVals = reshape(uVals.',[2 1 size(uVals,1)]);
        end


        function Re = elementalRotationMatrix(obj,deformedCoords)
                midPoints= 0.5*[deformedCoords(1,1)+deformedCoords(4,1),deformedCoords(1,2)+deformedCoords(4,2);
                             deformedCoords(2,1)+deformedCoords(3,1),deformedCoords(2,2)+deformedCoords(3,2)];
                m = [midPoints(2,1)-midPoints(1,1),midPoints(2,2)-midPoints(1,2)];
                mx = m(1); my = m(2);   
                Re = [mx, -my; my, mx] / sqrt(mx^2 + my^2);
        end
    end
end
