classdef CohesiveSeparationComputer < handle
    properties (Access = public)
        cohesiveMesh
        subMesh
        globalSeparations % nCohElem x nMidNodesPerElem(2) x nSeparations(2)

    end

    methods (Access = public)

        function obj = CohesiveSeparationComputer(cParams)
            obj.init(cParams);
            obj.createSubMesh();
            obj.ComputeGlobalSeparations()
        end
        










        function globalSeps = ComputeGlobalSeparations(obj)
        
            connec = obj.subMesh.connec;   % (nElem × 4)
            coords = obj.subMesh.coord;    % (nNodes × ndim)
        
            B = [-1 0 0 1;
                  0 -1 1 0];                         % (2×4)
        
            
        end




   end


    methods (Access = private)
        function init(obj,cParams)
            obj.cohesiveMesh = cParams.cohesiveMesh;
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


    end







end
