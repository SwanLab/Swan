classdef CohesiveMesh < handle
    
    properties (Access = public)
        mesh
        isCohesive

        newCoord
        newConnec

        nnodeCohesive
        pairsMatrix
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
            
        end
        
        function baseMesh(obj)
            obj.mesh = UnitQuadMesh(3,3);
        end

        function duplicator(obj)


            ymin = min(obj.mesh.coord(:,2));
            obj.isCohesive = abs(obj.mesh.coord(:,2)) == ymin; 
            obj.nnodeCohesive = sum(obj.isCohesive == 1) * 2;
            
            obj.newCoord = obj.mesh.coord;


            duplicated = obj.newCoord(obj.isCohesive == 1, :) + [0 , 0.1];
            obj.newCoord = [obj.newCoord; duplicated];

            obj.pairsMatrix = [find(obj.isCohesive == 1) , (size(obj.isCohesive,1)+1:1:obj.nnodeCohesive/2+size(obj.isCohesive,1))'];

            obj.isCohesive = [obj.isCohesive; ones(obj.nnodeCohesive/2,1)];

                scatter(obj.newCoord(:,1), obj.newCoord(:,2), 'filled');

        end

        function updateConnec(obj)
            obj.newConnec = obj.mesh.connec;
            
            cohesiveConnec = [obj.pairsMatrix(1:end-1,1), obj.pairsMatrix(2:end,1), obj.pairsMatrix(2:end,2), obj.pairsMatrix(1:end-1,2)];

            for i = obj.pairsMatrix(:,1)
                
        
            end



        end

        
    end
    
end