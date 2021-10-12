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
        alphaRot
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
            obj.nx1    = 152;
            obj.nx2    = 152;            
            obj.nCells = 16;
        end
        
        function loadDataExperiment(obj)
            iteration = 216;
            s.fileName = [obj.fileName,num2str(iteration)];
            s.folderPath = fullfile(obj.filePath );
            w = WrapperMshResFiles(s);
            w.compute();
            w = obj.symmetrizeData(w);
            obj.m1M    = w.m1;
            obj.m2M    = w.m2;
            obj.qM     = w.q;
            obj.mesh   = w.mesh;
            obj.alphaM = w.alpha;%obj.interpolateAlpha(w.dataRes.AlphaGauss');
        end
        
       function [x1r,x2r] = rotateCoord(obj,x1,x2)            
            ca = cos(obj.alphaRot);
            sa = sin(obj.alphaRot);                      
            x1r = ca*x1 - sa*x2;
            x2r = sa*x1 + ca*x2;                        
        end
        
        function wC = symmetrizeData(obj,w)
            coordT = w.mesh.coord;
            coordT(:,1) = coordT(:,1) + 0;%3;            
            coordT(:,2) = coordT(:,2) + 0;%2;
         %   [x1,x2] = obj.rotateCoord(coordT(:,1),coordT(:,2));
    %        coordT(:,1) = x1;
    %        coordT(:,2) = x2;
            
            sR.coord = coordT;
            sR.connec  = w.mesh.connec;
            w.mesh = Mesh(sR);
            
            m = w.mesh;
            sM = obj.createMeshSymetrizer(m);
            wC.mesh = sM.computeSymmetricMesh();
            

            alpha = obj.interpolateAlpha(w.dataRes.AlphaGauss',m);
            
            wC.m1 = sM.symmetrizeScalarField(w.dataRes.DesignVar1);
            wC.m2 = sM.symmetrizeScalarField(w.dataRes.DesignVar2);
            wC.q  = sM.symmetrizeScalarField(w.dataRes.SuperEllipseExponent);
            wC.alpha = sM.symmetrizeVectorField(alpha);              
            
           % sQ.mesh = w.mesh;
           % sQ.symmetricLine.vector = [1;0];
            %sQ.symmetricLine.point = [0;2];            
            %sM2 = Symmetrizer(sQ);
            %sM2.compute();
            
%             wC.m1 = w.dataRes.DesignVar1;         
%             wC.m2 = w.dataRes.DesignVar2;         
%             wC.q  = w.dataRes.SuperEllipseExponent;  
%             wC.alpha = alpha;
%             wC.mesh = w.mesh;
            
        end
        
 
        
        function mS = createMeshSymetrizer(obj,mesh)
            s.mesh = mesh;
            s.symmetricLine.vector = [1;0];
            s.symmetricLine.point = [0,0];           
            mS = Symmetrizer(s);
        end
        
        function alphaP1 = interpolateAlpha(obj,alphaM,m)
            alphaP0 = alphaM';
            s.mesh        = m;
            s.femSettings.scale = 'MACRO';
            s.femSettings.mesh = m;
            s.quadratureOrder = 'LINEAR';
            s.fileName = [];
            filter = Filter_P1_Density(s);
            alphaP1(:,1) = filter.getP1fromP0(alphaP0(:,1));
            alphaP1(:,2) = filter.getP1fromP0(alphaP0(:,2));
        end
        
        function createBackgroundMesh(obj)
            x1 = linspace(min(obj.mesh.coord(:,1)),max(obj.mesh.coord(:,1)),obj.nx1);
            x2 = linspace(min(obj.mesh.coord(:,2)),max(obj.mesh.coord(:,2)),obj.nx2);
            
            
            FV.vertices = [obj.mesh.coord,zeros(size(obj.mesh.coord,1),1)];
            FV.faces    = obj.mesh.connec;
            FV2 = refinepatch(FV);
            FV2 = refinepatch(FV2);
            %FV2 = refinepatch(FV2);
         %   FV2 = refinepatch(FV2);

            
            s.coord = FV2.vertices(:,1:2);
            s.connec = FV2.faces;
            m = Mesh(s);
            
            obj.backgroundMesh = m;
            
            %x1T = repmat(x1,obj.nx2,1);
            %x2T = repmat(x2',1,obj.nx1);
            %[x1T,x2T] = meshgrid(x1,x2); 

%             x1T = x1;
%             x2T = x2;
%             coord = [x1T(:),x2T(:)];            
%             xy = coord;
%             connec   = delaunay(xy(:,1),xy(:,2));
%             s.coord  = [xy(:,1),xy(:,2)];
%             s.connec = connec;
%             obj.backgroundMesh = Mesh(s);
            
            
%             x1min = min(x1);
%             x1max = max(x1);
%             x2min = min(x2);
%             x2max = max(x2);
%             [coordinates, nodes,nel,nnode] = MeshRectangular(x1max-x1min,x2max-x2min,obj.nx1,obj.nx2);
%             s.coord(:,1) = coordinates(:,1)+x1min;
%             s.coord(:,2) = coordinates(:,2)+x2min;
%             s.connec = nodes;
%             obj.backgroundMesh = Mesh(s);   
            
%              [xv,yv] = meshgrid(x1,x2); 
%              [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');             
%              s.coord  = V(:,1:2);
%              s.connec = F;
%              obj.backgroundMesh = Mesh(s);            
        end
        
        function createOrientation(obj)
            x1 = (obj.alphaM(:,1));
            x2 = (obj.alphaM(:,2));
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
            s.mesh               = obj.mesh;
            d = Dehomogenizer(s);
            d.compute();
            d.plot();
        end
        
    end
    
end