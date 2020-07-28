classdef PlottingToyUnfittedExample < testNotShowingError
    
    properties (SetAccess = protected, GetAccess = private)
        coord
        connec
        levelSet
    end
    
    properties (Access = private)
       backgroundMesh 
       boundaryMeshes
       unfittedMesh
    end
    

    methods (Access = protected)
        
        function compute(obj)
            obj.computeBackgroundMesh();
            obj.computeBoundaryMeshes();
            obj.computeUnfittedMesh();
            obj.plotUnfittedMesh();            
        end
        
        function hasPassed = hasPassed(obj)
            d = load(obj.testName);     
            itIs = isequal(obj.unfittedMesh,d.unfittedMesh);            
            hasPassed = itIs;
        end        
        
    end    
    
    methods (Access = private)
       
       function computeBackgroundMesh(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
            obj.backgroundMesh = Mesh(s);
        end
        
        function computeBoundaryMeshes(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.dimension = 1:obj.backgroundMesh.ndim;
            bC = BoundaryMeshCreatorFromRectangularBox(s);
            obj.boundaryMeshes = bC.create();                        
        end
        
        function computeUnfittedMesh(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.boundaryMesh   = obj.boundaryMeshes;
            uM = UnfittedMesh(s);
            uM.compute(obj.levelSet);
            obj.unfittedMesh = uM;
        end
        
        function plotUnfittedMesh(obj)
            obj.unfittedMesh.plot();
        end
    
    end
    
end