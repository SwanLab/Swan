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
            
            bMesh = obj.unfittedMesh.boundaryCutMesh.mesh;
            
            quad = Quadrature.set(bMesh.type);
            quad.computeQuadrature('QUADRATIC');
            
            
            xCutB = bMesh.computeXgauss(quad.posgp);
            p = plot(xCutB(1,:),xCutB(2,:),'s','MarkerSize',5); 
            color = [0.9290    0.6940    0.1250];
            p.MarkerEdgeColor = color;
            p.Color = color;
            p.MarkerFaceColor = color;
            
      
            bMesh = obj.unfittedMesh.innerCutMesh.mesh;   
            quad = Quadrature.set(bMesh.type);
            quad.computeQuadrature('QUADRATIC');                  
            xCutB = bMesh.computeXgauss(quad.posgp);
            p = plot(xCutB(1,:),xCutB(2,:),'s','MarkerSize',5);  
            color = [0.4940    0.1840    0.5560];
            p.MarkerEdgeColor = color;
            p.Color = color;
            p.MarkerFaceColor = color;            
            
            bMesh = obj.unfittedMesh.innerMesh.mesh;   
            quad = Quadrature.set(bMesh.type);
            quad.computeQuadrature('LINEAR');                  
            xCutB = bMesh.computeXgauss(quad.posgp);
            color = [0.8500    0.3250    0.0980];
            p = plot(xCutB(1,:),xCutB(2,:),'s','MarkerSize',5);  
            p.MarkerEdgeColor = color;
            p.Color = color;
            p.MarkerFaceColor = color;         
            
            bMesh = obj.unfittedMesh.unfittedBoundaryMesh;
            meshes = bMesh.getActiveMesh();
            for iMesh = 1:numel(meshes)
                mesh = meshes{iMesh};
                if ~isempty(mesh.innerMesh)
                bMesh = mesh.innerMesh.mesh;    
                quad = Quadrature.set(bMesh.type);
                 quad.computeQuadrature('QUADRATIC');                        
                xCutB = bMesh.computeXgauss(quad.posgp);
                hold on
                p = plot(xCutB(1,:),xCutB(2,:),'s','MarkerSize',5); 
                color = [0.9290    0.6940    0.1250];
                p.MarkerEdgeColor = color;
                p.Color = color;
                p.MarkerFaceColor = color;                
                end
                
                if ~isempty(mesh.innerCutMesh)
                bMesh = mesh.innerCutMesh.mesh;        
                quad = Quadrature.set(bMesh.type);
                quad.computeQuadrature('QUADRATIC');                  
                xCutB = bMesh.computeXgauss(quad.posgp);
                hold on
                p = plot(xCutB(1,:),xCutB(2,:),'s','MarkerSize',5); 
                color = [0.9290    0.6940    0.1250];
                p.MarkerEdgeColor = color;
                p.Color = color;
                p.MarkerFaceColor = color;                
                end
                
            end
            
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