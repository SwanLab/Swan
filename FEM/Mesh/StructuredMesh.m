classdef StructuredMesh < Mesh
    
    properties (Access = public)
        nx
        ny
        x
        y        
    end
    
    properties (Access = private)
       xv
       yv
    end
   
    methods (Access = public)
        
        function obj = StructuredMesh(cParams)
            obj.init(cParams);
            s.coord  = obj.createCoordinates();
            s.connec = obj.createConnectivities();
            obj.create(s);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)            
            [obj.x,obj.y] = meshgrid(cParams.x,cParams.y);
            obj.nx = length(obj.x);
            obj.ny = length(obj.y);
        end
        
        function coord = createCoordinates(obj)
            coord(:,1) = reshape(obj.x',1,[]);
            coord(:,2) = reshape(obj.y',1,[]);
        end
                
        function connec = createConnectivities(obj)
            Nx = obj.nx;
            Ny = obj.ny;
            sw = obj.createVertex(1:Nx-1,1:Ny-1);
            se = obj.createVertex(2:Nx,1:Ny-1);
            ne = obj.createVertex(2:Nx,2:Ny);
            nw = obj.createVertex(1:Nx-1,2:Ny);
            connec = [sw se ne nw];
        end        
        
        function v = createVertex(obj,xv,yv)
            Ny = obj.ny;
            v = bsxfun(@(x,y) x + Ny*(y-1),xv',yv);
            v = v(:);
        end        
        
    end
    
end