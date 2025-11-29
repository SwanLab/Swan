classdef CohesiveMesh < handle
    
    properties (Access = public)
        mesh
        isCohesive

        newCoord
        newConnec

        listNodeCohesive
        listElemCohesive

        nnodeCohesive
        pairsMatrix

        separation
    end
    
    properties (Access = private)

    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CohesiveMesh()
            obj.init()
            obj.baseMesh()
            obj.duplicator()
            obj.updateConnec()        
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.separation = 0.1;
        end
        
        function baseMesh(obj)
            obj.mesh = UnitQuadMesh(2,2);
        end

        function duplicator(obj)

            ymin = min(obj.mesh.coord(:,2));
            obj.isCohesive = abs(obj.mesh.coord(:,2)) == ymin; 
            obj.nnodeCohesive = sum(obj.isCohesive == 1) * 2;
            obj.listNodeCohesive = find(obj.isCohesive == 1);

            obj.newCoord = obj.mesh.coord;

            dispV      = [0, 1];
            duplicated = obj.newCoord(obj.isCohesive == 1, :) + obj.separation*dispV;
            obj.newCoord = [obj.newCoord; duplicated];

            obj.pairsMatrix = [obj.listNodeCohesive , (size(obj.isCohesive,1)+1:1:obj.nnodeCohesive/2+size(obj.isCohesive,1))'];
            obj.pairsMatrix = [obj.pairsMatrix , dispV(1)*ones(obj.nnodeCohesive/2,1), dispV(2)*ones(obj.nnodeCohesive/2,1)];

            obj.isCohesive = [obj.isCohesive; ones(obj.nnodeCohesive/2,1)];

                scatter(obj.newCoord(:,1), obj.newCoord(:,2), 'filled');

        end

        function updateConnec(obj)
            obj.newConnec = obj.mesh.connec;
            
            isCohesiveElem = any(obj.isCohesive(obj.newConnec),2);
            obj.listElemCohesive = find(isCohesiveElem == 1);

            % No se com identificar quins s'han d'intercanviar d'una forma
            % eficient

            cohesiveConnec = [obj.pairsMatrix(1:end-1,1), obj.pairsMatrix(2:end,1), obj.pairsMatrix(2:end,2), obj.pairsMatrix(1:end-1,2)];
            obj.newConnec = [obj.newConnec; cohesiveConnec];

        end
        
    end
    
end