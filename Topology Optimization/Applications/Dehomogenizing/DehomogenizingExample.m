classdef DehomogenizingExample < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        m1M
        m2M
        qM
        alphaM
        mesh
        backgroundMesh
        theta
        superEllipse
        cellLevelSetParams
    end
    
    properties (Access = private)
        filePath
        fileName
        nCells
        nx1
        nx2
    end
    
    methods (Access = public)
        
        function obj = DehomogenizingExample()
            obj.init();
            obj.loadDataExperiment();
            obj.createBackgroundMesh();
            obj.createSuperEllipseParams();
            obj.createOrientation();
            obj.createLevelSetCellParams();
            obj.dehomogenize();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.filePath = '/home/alex/git-repos/Swan/Topology Optimization/Applications/Dehomogenizing/Example/';
            obj.fileName = 'CantileverSymmetricWithoutFixing';
            obj.nx1    = 352;
            obj.nx2    = 352;            
            obj.nCells = 55;
        end
        
        function loadDataExperiment(obj)
            iteration = 216;
            s.fileName = [obj.fileName,num2str(iteration)];
            s.folderPath = fullfile(obj.filePath );
            w = WrapperMshResFiles(s);
            w.compute();
            obj.m1M    = w.dataRes.DesignVar1;
            obj.m2M    = w.dataRes.DesignVar2;
            obj.qM     = w.dataRes.SuperEllipseExponent;
            obj.mesh   = w.mesh;
            obj.alphaM = obj.interpolateAlpha(w.dataRes.AlphaGauss');
        end
        
        function alphaP1 = interpolateAlpha(obj,alphaM)
            alphaP0 = alphaM';
            s.mesh        = obj.mesh;
            s.femSettings.scale = 'MACRO';
            s.femSettings.mesh = obj.mesh;
            s.quadratureOrder = 'LINEAR';
            s.fileName = [];
            filter = Filter_P1_Density(s);
            alphaP1(:,1) = filter.getP1fromP0(alphaP0(:,1));
            alphaP1(:,2) = filter.getP1fromP0(alphaP0(:,2));
        end
        
        function createBackgroundMesh(obj)
            x1 = linspace(min(obj.mesh.coord(:,1)),max(obj.mesh.coord(:,1)),obj.nx1);
            x2 = linspace(min(obj.mesh.coord(:,2)),max(obj.mesh.coord(:,2)),obj.nx2);
            x1T = repmat(x1,obj.nx2,1);
            x2T = repmat(x2',1,obj.nx1);
            coord = [x1T(:),x2T(:)];            
            xy = coord;
            x1 = xy(:,1);
            x2 = xy(:,2);
            connec   = delaunay(x1,x2);
            s.coord  = [x1,x2];
            s.connec = connec;
            obj.backgroundMesh = Mesh(s);
        end
        
        function createOrientation(obj)
            x1 = obj.interpolateFunction(obj.alphaM(:,1));
            x2 = obj.interpolateFunction(obj.alphaM(:,2));
            obj.theta = atan2(x2,x1);
        end
        
        function createSuperEllipseParams(obj)
            sE.m1 = obj.interpolateFunction(obj.m1M);
            sE.m2 = obj.interpolateFunction(obj.m2M);
            sE.q  = obj.interpolateFunction(obj.qM);
            obj.superEllipse = sE;
        end
        
        function vq = interpolateFunction(obj,v)
            X = obj.mesh.coord(:,1);
            Y = obj.mesh.coord(:,2);
            F = scatteredInterpolant(X,Y,v);
            xB = obj.backgroundMesh.coord(:,1);
            yB = obj.backgroundMesh.coord(:,2);
            vq = F(xB,yB);
        end        
        
        function createLevelSetCellParams(obj)
            s.type   = 'smoothRectangle';
            s.widthH = obj.superEllipse.m1;
            s.widthV = obj.superEllipse.m2;
            s.pnorm  = obj.superEllipse.q;
            s.ndim   = 2;
            obj.cellLevelSetParams = s;
        end
        
        function dehomogenize(obj)
            s.backgroundMesh     = obj.backgroundMesh;
            s.nCells             = obj.nCells;
            s.theta              = obj.theta;
            s.cellLevelSetParams = obj.cellLevelSetParams;
            d = Dehomogenizer(s);
            d.compute();
            d.plot();
        end
        
    end
    
end