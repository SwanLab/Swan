classdef TopOptTestTutorial < handle

    properties (Access = private)
        mesh
    end

    methods (Access = public)

        function obj = TopOptTestTutorial()
            obj.init()
            obj.createMesh();
        end

    end

    methods (Access = private)

        function init(obj)

        end

        function createMesh(obj)
            x1      = linspace(0,1,50);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh(s);            
        end

        function createBoundaryConditions(obj)
            
        end

    end

end