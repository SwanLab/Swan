classdef Jump < handle

    properties (Access = private)
        jumpDim
        L
        jumpFun
    end

    properties (Access = private)
        uFun
        cohesiveMesh     
        ndimf 
    end

    methods (Access = public)

        function obj = Jump(cParams)     
            obj.init(cParams);
            obj.createJumpFunction();
            obj.computeGlobalSeparationMatrix();
            obj.updateJumpValues(obj.uFun); 
        end
        function updateJumpValues(obj,uIn)          
            % uInVec    = reshape(uIn.fValues',[uIn.nDofs 1]); ? 
            R = obj.computeRotationMatrix; % 2x2xnElem
            connJump  = CohesiveMesh.lineMesh.connec;
            nNodeJump = obj.cohesiveMesh.lineMesh.nnode;
            fValues   = zeros(nNodeJump,obj.jumpDim); % nNode x 2
            for i = 1:uIn.ndimf
                uVals = obj.computeDispDofs(uIn, i);
                jumpi = obj.L(:,:,i) * R * uVals; %2x1xnElem
                localNode = 1 + ~(i<=2 || i>=7);
                globalNode = connJump(:,localNode);
                fValues(globalNode,:) = fValues(idx,:) + jumpi;
            end
            obj.jumpFun.setFValues(fValues);
        end
        function fV = evaluate(obj,xV)
            fV = obj.jumpFun.evaluate(xV);
        end
        
        
        function N = computeShapeFunctions(obj,xV) % pel test
            R = obj.computeRotationMatrix(obj.uFun); % 2x2xnElem
            L = obj.L;          % 2x2x8
            Njump = obj.jumpFun.computeShapeFuncion(xV); 
            
            N
            
            %N = N*L*R;
            % 2 x 8 x nElem resultat?
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
                    % % L -- ndofJump x ndofu
                    % nCohElem     = length(obj.cohesiveMesh.listCohesiveElems);
                    % nJumpPerElem = obj.jumpDim * obj.cohesiveMesh.lineMesh.nnodeElem;
                    % nDofU        = obj.cohesiveMesh.mesh.nnodes * obj.cohesiveMesh.mesh.ndim;
                    % nDofJump     = nJumpPerElem* nCohElem;
                    % ndim         = obj.cohesiveMesh.mesh.ndim;
                    % L_full       = zeros(nDofJump, nDofU);
                    % 
                    % Le = [-1,0,0,0,0,0,1,0;
                    %        0,-1,0,0,0,0,0,1;
                    %        0,0,-1,0,1,0,0,0;
                    %        0,0,0,-1,0,1,0,0];
                    % 
                    % for j = 1:nCohElem
                    %     elem   = obj.cohesiveMesh.listCohesiveElems(j);
                    %     connec = obj.cohesiveMesh.mesh.connec(elem,:);
                    % 
                    %     dofsU    = reshape(((connec(:)-1)*ndim + (1:ndim)).', 1, []);
                    %     dofsJump = nJumpPerElem*(j-1)*ones(nJumpPerElem,1).'+(1:nJumpPerElem);
                    % 
                    %     L_full(dofsJump,dofsU) = Le;
                    % end
                    % obj.L = sparse(L_full);
            obj.L =cat(3, -repmat(eye(2),1,1,4), repmat(eye(2),1,1,4));
        end
        
        function Rall = computeRotationMatrix(obj,uIn) 
            nDofsU       = uIn.nDofs;
            nCohElem     = length(obj.cohesiveMesh.listCohesiveElems);
            nJumpPerElem = obj.jumpDim * obj.cohesiveMesh.lineMesh.nnodeElem;

            Rfull = zeros(nDofsU,nDofsU);
            Rall  = zeros(2,2,nCohElem);

            for j=1:nCohElem
                e = obj.cohesiveMesh.listCohesiveElems(j);
                connecMesh = obj.cohesiveMesh.mesh.connec(e,:);
                coordsMesh = obj.cohesiveMesh.mesh.coord(connecMesh',:);
                %dofsU = reshape(((connecMesh(:)-1)*obj.ndimf + (1:obj.ndimf)).', 1, []);
                disp  = uIn.fValues(connecMesh',:);
                    deformedCoords = coordsMesh + disp;
                Re = elementalRotationMatrix(deformedCoords);
                %ReBig              = kron(eye(nJumpPerElem), Re);
                Rall(:,:,j)        = Re;
                %Rfull(dofsU,dofsU) = ReBig;
            end
            %R = sparse(Rfull);
        end
        
        function uVals = computeDispDofs(uIn,i) % 2x1xnElem, Verificar!
            fValuesU = uIn.fValues;
            conn = uIn.mesh.connec;
            dof = [uIn.ndimf*conn-1  uIn.ndimf*conn];
            dof = reshape(dof.',uIn.nDofsElem,[]).'; 
            dof = dof(:,i);
            uVals = fValuesU(dof,:);
        end
    
        function Re = elementalRotationMatrix(deformedCoords)
                midPoints= 0.5*[deformedCoords(1,1)+deformedCoords(4,1),deformedCoords(1,2)+deformedCoords(4,2);
                             deformedCoords(2,1)+deformedCoords(3,1),deformedCoords(2,2)+deformedCoords(3,2)];
                m = [midPoints(2,1)-midPoints(1,1),midPoints(2,2)-midPoints(1,2)];
                mx = m(1); my = m(2);   
                Re = [mx, -my; my, mx] / sqrt(mx^2 + my^2);
        end
    end
end
