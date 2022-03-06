classdef PerformanceTest < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        dim
        type
        scale
        mesh
        length, height
        dirichlet
        neumann
    end
    
    methods (Access = public)
        
        function obj = PerformanceTest(cParams)
            obj.init(cParams)
        end

        function sol = compute(obj)
            obj.createCantileverBeam();
            obj.computeBoundaryConditions();
            sol = obj.computeSolution();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.dim    = cParams.dim;
            obj.length = cParams.length;
            obj.height = cParams.height;
        end
        
        function createCantileverBeam(obj)
            s.dim    = obj.dim;
            s.len    = obj.length;
            s.height = obj.height;
            beam     = CantileverBeam(s);
            obj.mesh = beam.create(5,5);
        end
        
        function computeBoundaryConditions(obj)
            coords = obj.mesh.coord;
            npnod = size(coords,1);
            dirichNodes = find(coords(:,1) == 0 );
            ndirich = size(dirichNodes,1);
            dir = [dirichNodes,   ones(ndirich,1), zeros(ndirich,1);
                   dirichNodes, 2*ones(ndirich,1), zeros(ndirich,1)];
            neu   = [npnod, 2, 0.0001];
            obj.dirichlet = dir;
            obj.neumann   = neu;
        end
        
        function fem = computeSolution(obj)
            p = obj.createFEMparameters();
            fem = NewFEM.create(p);
            fem.solve();
        end
        
        function p = createFEMparameters(obj)
            p.dim       = obj.dim;
            p.type      = 'ELASTIC';
            p.scale     = 'MACRO';
            p.mesh      = obj.mesh;
            p.dirichlet = obj.dirichlet;
            p.pointload = obj.neumann;
        end
    end
    
end