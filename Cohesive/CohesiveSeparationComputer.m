classdef CohesiveSeparationComputer < handle
    properties (Access = public)
        cohesiveMesh
        subMesh
        globalSeparations % nCohElem x nMidNodesPerElem(2) x nSeparations(2)
        lagrangianSeparation
        u
    end

    methods (Access = public)

        function obj = CohesiveSeparationComputer(cParams)
            obj.init(cParams);
            obj.createSubMesh();
            obj.globalSeparations = obj.ComputeGlobalSeparations();
            obj.lagrangianSeparation = obj.ComputeLagrangianSeparation();
        end
        
        function globalSeps = ComputeGlobalSeparations(obj)
        
            connec = obj.cohesiveMesh.mesh.connec;   % (nElem × 4)
            coords = obj.cohesiveMesh.mesh.coord;    % (nNodes × ndim)
        
            L = [-1 0 0 1;
                  0 -1 1 0];                         % (2×4)
        
            nelem = length(obj.cohesiveMesh.listCohesiveElems);
            globalSeps = zeros(nelem,2,2);

            for i = 1:nelem
                e = obj.cohesiveMesh.listCohesiveElems(i);
                R = obj.rotationMatrix(i);
                Xe = coords(connec(e,:)',:);
                Ue = obj.u.fValues(connec(e,:)',:);
                globalSeps(i,:,:) = L*(Xe+Ue)*R';
            end
            
        end

   end


    methods (Access = private)
        function init(obj,cParams)
            obj.cohesiveMesh = cParams.cohesiveMesh;
            obj.u = cParams.u;
        end

        function createSubMesh(obj)

            coord = obj.cohesiveMesh.mesh.coord;
            nMidNodes = size(obj.cohesiveMesh.pairsMatrix,1);
            
            midCoord= (coord(obj.cohesiveMesh.listNodeCohesive,:)+ ...
                coord(obj.cohesiveMesh.pairsMatrix(:,2),:))/2;
            
                % midCoord = midCoord(:,1);
            
            midConnec =  [(1:nMidNodes-1)' (2:nMidNodes)'];
            
                s.connec = midConnec;
                s.coord = midCoord;
                s.kFace = -1;
            obj.subMesh = Mesh.create(s);

        end

        function R = rotationMatrix(obj,i)
            connec = obj.subMesh.connec(i,:);
            coords = obj.subMesh.coord(connec',:);
            
            m = [coords(2,1)-coords(1,1),coords(2,2)-coords(1,2)];
                mx = m(1);
                my = m(2);
            R = [mx, -my; my, mx] / sqrt(mx^2 + my^2);

        end



        function lagrangianSeparation = ComputeLagrangianSeparation(obj)
            




            %                   AMB MITJA

            % globalSeps = reshape(obj.globalSeparations, [], 2)';
            % tmp = globalSeps(:,2:(size(globalSeps,2)-1));
            % pairs = tmp(:,2:2:end);
            % odds = tmp(:,1:2:end);
            % tmp = (pairs+odds)/2;
            % meanSeps = [globalSeps(:,1), tmp, globalSeps(:,end)]';
            % fValues = meanSeps;
            % lagrangianSeparation = LagrangianFunction.create(obj.subMesh,2,'P1');
            % lagrangianSeparation.setFValues(fValues);





            %                   AMB FUNCIO DISCONTINUA

            fValues = globalSeps';
            lagrangianSeparation = LagrangianFunction.create(obj.subMesh,2,'P1D');
            lagrangianSeparation.setFValues(fValues);


        end










    end


end
