classdef HarmonicProjectionExample < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        filePath
        fileName
        iteration 
        experimentData
        mesh
        orientationAngle
        orientationAngleGauss
    end
    
    methods (Access = public)
        
        function obj = HarmonicProjectionExample()
            obj.init();
            obj.loadDataExperiment();
            obj.storeMesh();
            obj.storeOrientationAngle();
            obj.plotOrientation();
            obj.project()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
           obj.filePath = '/home/alex/git-repos/Swan/Topology Optimization/Applications/Dehomogenizing/ExampleLShape/';
           obj.fileName = 'LshapeCoarseSuperEllipseDesignVariable';
           obj.iteration = 665;           
        end
        
        function loadDataExperiment(obj)            
            s.fileName = [obj.fileName,num2str(obj.iteration)];
            s.folderPath = fullfile(obj.filePath );
            w = WrapperMshResFiles(s);
            w.compute();
            obj.experimentData = w;
        end        
        
        function storeMesh(obj)
            d = obj.experimentData;
            obj.mesh = d.mesh;
        end
        
        function storeOrientationAngle(obj)
            d = obj.experimentData;
            alpha  = d.dataRes.AlphaGauss;
            alpha = obj.interpolateOrientationAngle(alpha);
            obj.orientationAngle = alpha;
        end
        
        function alphaP1 = interpolateOrientationAngle(obj,alphaP0)
            s.mesh        = obj.mesh;
            s.femSettings.scale = 'MACRO';
            s.femSettings.mesh = obj.mesh;
            s.quadratureOrder = 'LINEAR';
            s.fileName = [];
            filter = Filter_P1_Density(s);
            alphaP1(:,1) = filter.getP1fromP0(alphaP0(:,1));
            alphaP1(:,2) = filter.getP1fromP0(alphaP0(:,2));
        end     

        function plotOrientation(obj)
             figure()
             x = obj.mesh.coord(:,1);
             y = obj.mesh.coord(:,2);
             t  = obj.orientationAngle;
             ct = t(:,1);
             st = t(:,2);
             quiver(x,y,ct,st)
        end     

        function project(obj)
            s.mesh = obj.mesh;
            s.orientationAngle = obj.orientationAngle;
            h = HarmonicProjection(s);
            h.project();
        end
        
    end
    
end