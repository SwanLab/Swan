classdef OptimalSuperEllipsePrinter < handle
    
    properties (Access = private)
        mesh
        meshBackground
        levelSet
        mx
        my     
        iter
    end
    
    properties (Access = private)
        mxV
        myV        
        fileName
    end
    
    methods (Access = public)
        
        function obj = OptimalSuperEllipsePrinter()
            obj.init();
            obj.compute();
        end  
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.mxV = [0.8462 0.99 0.01 0.99 0.30108 0.1 0.2 0.2  0.85];
            obj.myV = [0.8462 0.2  0.01 0.99 0.72687 0.1 0.6 0.95 0.85];
            obj.fileName = 'OptimaSuperEllipseMicroEllipse';
          %  iter = 0;
        end
        
        function compute(obj)            
            for iTest = 1:length(obj.mxV)
                obj.iter = iTest;             
                
                rho = SuperEllipseParamsRelator.rho(obj.mxV(obj.iter),obj.myV(obj.iter),32)
                
                obj.createMeshBackground();
                obj.createLevelSet();
                obj.createMesh();
                obj.print();
            end
        end
        
        function createMeshBackground(obj)
            run('RVE_Square_Triangle_FineFine')
            a.connec = connec(:, 2:end);
            a.coord  = coord(:, 2:3);
            m = Mesh.create(a);
            %obj.testName = 'RVE_Square_Triangle_Fine';
            obj.meshBackground = m; 
        end
        
        function createLevelSet(obj)
            sM.coord  = obj.meshBackground.coord;
            sM.connec = obj.meshBackground.connec;
            s.mesh = Mesh_Total(sM);
            s.widthH = obj.mxV(obj.iter);
            s.widthV = obj.myV(obj.iter);
            %s.pnorm  = obj.computeSmoothingExponent();
            s.pnorm  = 2;
            s.type = 'smoothRectangle';
            s.levelSetCreatorSettings = s;
            s.type = 'LevelSet';
            s.scalarProductSettings.epsilon = 1;
            obj.levelSet = LevelSet(s);
        end
        
        function q = computeSmoothingExponent(obj)
            s.m1 = obj.mxV(obj.iter);
            s.m2 = obj.myV(obj.iter);
            s.type = 'Optimal';
            qComputer = SmoothingExponentComputer.create(s);
            q = qComputer.compute();
        end
        
        function createMesh(obj)
            s.fileName = obj.fileName;
            s.levelSet = obj.levelSet.value;
            s.meshBackground = obj.meshBackground;
            mCreator = MeshCreatorFromLevelSetWithMMG(s);
            obj.mesh = mCreator.create();
        end
        
        function print(obj)
            outputName = [obj.fileName,'Print',num2str(obj.iter)];
            s.mesh       = obj.mesh;
            s.outPutName = outputName;
            printer = SuperEllipsePrinter(s);
            printer.print();
            printer.captureImage();
        end
        
    end
    
end