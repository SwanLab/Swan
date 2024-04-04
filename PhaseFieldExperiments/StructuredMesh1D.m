classdef StructuredMesh1D < handle
    
    properties (Access = public)
        nx
        x
        mesh
    end
    
    properties (Access = private)
       xv
    end
   
    methods (Access = public)
        
        function obj = StructuredMesh1D(cParams)
            obj.init(cParams);
            s.coord  = obj.createCoordinates();
            s.connec = obj.createConnectivities();
            obj.createMesh(s);
        end

        function [xL,cells] = obtainLocalFromGlobalCoord(obj,xG)
            s.mesh     = obj;
            s.points.x = xG;
            cellFinder = CellFinderInStructuredMesh1D(s);
            xL    = cellFinder.naturalCoord;
            cells = cellFinder.cells;            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            [obj.x] = ndgrid(cParams.x);
            obj.nx = length(obj.x);

        end
        
        function coord = createCoordinates(obj)
            coord(:,1) = reshape(obj.x',1,[]);
        end
        
        function connec = createConnectivities(obj)
            Nx = obj.nx;
            e = (2:Nx)';
            w = (1:Nx-1)';
            connec = [w e];
        end

        function createMesh(obj,s)
            s.kFace = 0;
            obj.mesh = Mesh.create(s);
        end
        
    end
    
end