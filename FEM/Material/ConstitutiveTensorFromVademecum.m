classdef ConstitutiveTensorFromVademecum < handle
    
    properties (Access = private)
        fileName
        loadVariables
        domain
        mesh
        CtensorValues
    end
    
    methods (Access = public)
        
        function obj = ConstitutiveTensorFromVademecum(cParams)
            obj.init(cParams);
            obj.loadVademecumVariables();
            obj.obtainVademecumDomain();
            obj.createMesh();
            obj.obtainConstitutiveTensorValues()
        end
        
        function C = computeCtensor(obj,x)
            nstre = size(obj.CtensorValues,1);
            xG = obj.mesh.xG;
            yG = obj.mesh.yG;
            C = zeros(nstre,nstre,size(x,1));
            for i = 1:nstre
                for j = 1:nstre
                    cv = squeeze(obj.CtensorValues(i,j,:,:));
                    c  = interp2(xG,yG,cv,x(:,1),x(:,2));
                    c2 = obj.interpolate(xG,yG,cv,x(:,1),x(:,2));                    
                    C(i,j,:) = c;
                end
            end
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
        
        function obtainConstitutiveTensorValues(obj)
            var = obj.loadVariables.variables;
            for imx = 1:length(obj.domain.mxV)
                for imy = 1:length(obj.domain.myV)
                    C(:,:,imx,imy) = var{imx,imy}.Ctensor;
                end
            end
            obj.CtensorValues = C;            
        end
        
        function z = interpolate(obj,xG,yG,zG,x,y)
            s.xG = xG;
            s.yG = yG;
            s.zG = zG;
            int = Interpolator(s);
            z = int.interpolate(x,y);
        end
        
    end
    
    
end
