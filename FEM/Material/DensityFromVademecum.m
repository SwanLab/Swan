classdef DensityFromVademecum < handle
    
    properties (Access = private)
        fileName
        loadVariables
        domain
        mesh
        DensityValues
    end
    
    methods (Access = public)
        
        function obj = DensityFromVademecum(cParams)
            obj.init(cParams);
            obj.loadVademecumVariables();
            obj.obtainVademecumDomain();
            obj.createMesh();
            obj.obtainDensityValues()
        end
        
        function rho = computeDensity(obj,x)
            xG = obj.mesh.xG;
            yG = obj.mesh.yG;
            rho = zeros(size(x,1),1);
            d = squeeze(obj.DensityValues);
            rho(:,1) = interp2(xG,yG,d,x(:,1),x(:,2));
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
        end
        
        function loadVademecumVariables(obj)
            matFile   = [obj.fileName,'.mat'];
            file2load = fullfile('Output',obj.fileName,matFile);
            v = load(file2load);
            obj.loadVariables = v.d;
        end
        
        function obtainVademecumDomain(obj)
            obj.domain.mxV = obj.loadVariables.domVariables.mxV;
            obj.domain.myV = obj.loadVariables.domVariables.myV;
        end
        
        function createMesh(obj)
            obj.createMeshGrid()
            obj.createCoordinates();
            obj.createConnectivities();
        end
        
        function createMeshGrid(obj)
            [xG,yG] = meshgrid(obj.domain.mxV,obj.domain.myV);
            obj.mesh.xG = xG;
            obj.mesh.yG = yG;
        end
        
        function createCoordinates(obj)
            x = reshape(obj.mesh.xG',1,[])';
            y = reshape(obj.mesh.yG',1,[])';
            obj.mesh.coord(:,1) = x;
            obj.mesh.coord(:,2) = y;
        end
        
        function createConnectivities(obj)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            obj.mesh.connec = delaunay(x,y);
        end
        
        function obtainDensityValues(obj)
            var = obj.loadVariables.variables;
            for imx = 1:length(obj.domain.mxV)
                for imy = 1:length(obj.domain.myV)
                    rho(:,:,imx,imy) = var{imx,imy}.volume;
                end
            end
            obj.DensityValues = rho;
        end
        
    end
    
    
end
