classdef testingInterpDiscCont < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)

    end
    
    properties (Access = private)
        experimentData
        filePath
        fileName
        iteration
        mesh
        orientationAngle
    end
    
    methods (Access = public)
        
        function obj = testingInterpDiscCont()
            obj.init()
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
%         function init(obj)
%             obj.filePath = '/home/alex/git-repos/Swan/Topology Optimization/Applications/Dehomogenizing/ExampleLShape/';
%             obj.fileName = 'LshapeCoarseSuperEllipseDesignVariable';
%             obj.iteration = 665;
%             obj.loadDataExperiment();
%         end
% 
%         function loadDataExperiment(obj)
%             s.fileName = [obj.fileName,num2str(obj.iteration)];
%             s.folderPath = fullfile(obj.filePath );
%             w = WrapperMshResFiles(s);
%             w.compute();
%             obj.experimentData = w;
%         end
% 
%         function createMesh(obj)
%             d = obj.experimentData;
%             obj.mesh = d.mesh;
%         end
% 
%         function computeOrientationAngle(obj)
%             d = obj.experimentData;
%             alpha0  = d.dataRes.AlphaGauss;
%             alpha(:,1) = obj.interpolateOrientationAngle(alpha0(:,1));
%             alpha(:,2) = obj.interpolateOrientationAngle(alpha0(:,2));
% 
% 
%             theta(:,1) = atan2(alpha(:,1),alpha(:,2));  
% 
%             %obj.plotOrientation(theta,1);
%             alpha = obj.projectInUnitBall(alpha);
%             theta(:,1) = atan2(alpha(:,1),alpha(:,2));
%             %obj.plotOrientation(theta,1);
%             obj.orientationAngle = theta;
%        end   


        function init(obj)
        end        

        function createMesh(obj)
            connec = [1 2 4;
                1 4 3];
            coord = [0 0;
                1 0;
                0 1;
                1 1];
            sC.connec = connec;
            sC.coord  = coord;
            mC = Mesh(sC);       
            obj.mesh = mC;
        end
% 
        function computeOrientationAngle(obj)
            obj.orientationAngle  = pi/180*[200;170;60;0];
        end
    
        function compute(obj)

            obj.createMesh();
            obj.computeOrientationAngle();
            beta = obj.orientationAngle;
            mC = obj.mesh;
            tC    = obj.createUnitOrientedVector(beta);

            close all
            
            
            figure()
            hold on
            mC.plot()
            obj.plotOrientation(mC,tC,'b')
            
            alpha = beta/2;
            tC = obj.createUnitOrientedVector(alpha);
            obj.plotOrientation(mC,tC,'r')
            
            mD = mC.createDiscontinousMesh();

            
            s.connec = mC.connec;
            s.type   = mC.type;
            s.fNodes = tC;
            tcF = FeFunction(s);
            tD = tcF.computeDiscontinousField();

       
            
            figure()
            hold on
            mD.plot()
            obj.plotOrientation(mD,tD,'r')
            
            s.meshCont   = mC;
            s.meshDisc   = mD;
            s.fieldDisc  = tD;
            
            sC = SymmetricContMapCondition(s);
            c  = sC.computeCondition();
            
            
            
        end

        function createDiscontinousMesh(obj)

        end

        function vI = interpolateOrientationAngle(obj,v0)
            s.mesh    = obj.mesh;
            s.fValues = v0;
            p = PieceWiseConstantFunction(s);
            vI = p.projectToLinearNodalFunction();
        end      

        function vP = projectInUnitBall(obj,v)
            u = UnitBallProjector([]);
            vP = u.project(v);
        end        
        
    end
    
    methods (Access = private, Static)
        
        function v = createUnitOrientedVector(alpha)
            v = [cos(alpha) sin(alpha)];
        end
        
        function plotOrientation(m,t,color)
            quiver(m.coord(:,1),m.coord(:,2),t(:,1),t(:,2),color)
        end
    end
    
end