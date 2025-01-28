classdef TestNaca < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        mesh
        levelSet
    end
    
    methods (Access = public)
        
        function obj = TestNaca()
            obj.createMesh();
            obj.createLevelSet();
            obj.createUnfitted();
        end
        
    end
    
    methods (Access = private)
                
        function createMesh(obj)      
            x1 = linspace(0,1,40);
            x2 = linspace(0,1,40);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.mesh = m;
        end
        
        function createLevelSet(obj)
            sLS.type        = 'Naca';
            sLS.xCoorCenter = 0;
            sLS.yCoorCenter = 0.5;
            sLS.radius      = 0.2;
            g               = GeometricalFunction(sLS);
            lsFun           = g.computeLevelSetFunction(obj.mesh);
            obj.levelSet    = lsFun.fValues;
        end
        
        function createUnfitted(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(obj.levelSet);       
            figure
        %    uMesh.plot()
            m = uMesh.createInnerMesh();
            m.plot();


            % to be improved
            points = m.coord;
            % Create an alphaShape for the convex inclusion
            inclusionShape = alphaShape(inclusionPoints, Inf); % Convex hull around inclusion

            % Perform Delaunay triangulation
            dt = delaunayTriangulation(points);

            triangles = dt;

            % Compute the centroids of each tetrahedron
            centroids = (points(triangles(:,1), :) + ...
                points(triangles(:,2), :) + ...
                points(triangles(:,3), :) )/ 3;

             % Determine which triangles are outside the inclusion
             insideInclusion = inShape(inclusionShape, centroids(:, 1), centroids(:, 2));

             % Keep only triangles outside the inclusion
             trianglesFiltered = triangles(~insideInclusion, :);

             % Plot results
             figure;
             hold on;

             % Plot the Delaunay triangulation (filtered)
             figure;
             triplot(trianglesFiltered, points(:, 1), points(:, 2), 'b');


            % Delaunay Triangulation in 3D
            s.connec = delaunay(points);
            s.coord = m.coord;
            m1 = Mesh.create(s);
            figure
            plot(m1)
        end
        
    end
    
end