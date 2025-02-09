classdef PlottingSomeOptimalMicroStructures < handle
    
    properties (Access = private)
        mesh
        m1
        m2
        m1v
        m2v
        q
        levelSet
        postProcess
        outputName
        name
    end
    
    methods (Access = public)
        
        function obj = PlottingSomeOptimalMicroStructures()
            obj.init();
            obj.computePrinting();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.m1v = [0.1 0.2 0.2 0.85 0.99];
            obj.m2v = [0.1 0.6 0.95 0.85 0.99];
            obj.name = 'OptimalSuperEllipse';
        end
        
        function computePrinting(obj)
            for icase = 1:length(obj.m1v)
                obj.m1 = obj.m1v(icase);
                obj.m2 = obj.m2v(icase);
                obj.outputName = [obj.name,num2str(icase)];
                obj.q(icase) = obj.computeCornerSmoothingExponent();
                obj.createMesh();
                obj.createLevelSet(icase);
                obj.createUnfittedMesh();
                obj.createPostProcess();
                obj.print(icase);
            end
        end
        
        function createMesh(obj)
            fileName='RVE_Square_Triangle_FineFine';
            femReader = FemInputReader_GiD();
            reader    = femReader.read(fileName);
            s.coord   = reader.coord;
            s.connec  = reader.connec;
            obj.mesh  = Mesh.create(s);
        end
        
        function createLevelSet(obj,icase)
            s.type = 'LevelSet';
            s.mesh = obj.mesh;
            s.scalarProductSettings = obj.createScalarProductSettings;
            s.levelSetCreatorSettings = obj.createLevelSetCreatorSettings(icase);
            desVar = DesignVariable.create(s);
            obj.levelSet = desVar.value;
        end
        
        function createUnfittedMesh(obj)
            m = obj.mesh;
            s.unfittedType = 'INTERIOR';
            s.backgroundMesh = m.innerMeshOLD;
            s.boundaryMesh   = m.boxFaceMeshes;
            s.isInBoundary = false;
            uMesh = UnfittedMesh(s);
            uMesh.compute(obj.levelSet);
            uMesh.plot()
        end
        
        function qv = computeCornerSmoothingExponent(obj)
            s.m1 = obj.m1;
            s.m2 = obj.m2;
            s.type = 'Optimal';
            qComputer = SmoothingExponentComputer.create(s);
            qv = qComputer.compute();
        end
        
        function s = createLevelSetCreatorSettings(obj,icase)
            s.type = 'smoothRectangle';
            s.widthH = obj.m1;
            s.widthV = obj.m2;
            s.pnorm  = obj.q(icase);
        end
        
        function createPostProcess(obj)
            outName = obj.outputName;
            dI.mesh            = obj.mesh;
            dI.outName         = outName;
            dI.pdim            = '2D';
            dI.ptype           = 'MICRO';
            ps = PostProcessDataBaseCreator(dI);
            dB = ps.getValue();
            dB.printers = 'LevelSet';
            postCase    = 'LevelSet';
            obj.postProcess = Postprocess(postCase,dB);
        end
        
        function print(obj,isample)
            d.x = obj.levelSet;
            obj.postProcess.print(isample,d);
            i = isample;
            f = obj.outputName;
            outPutNameWithIter = [obj.outputName,num2str(i)];
            inputFileName = fullfile('Output',f,[f,num2str(i),'.flavia.res']);
            sI.fileName = f;
            sI.outPutImageName = outPutNameWithIter;
            sI.inputFileName = inputFileName;
            imageCapturer = GiDImageCapturer(sI);
            imageCapturer.capture();
        end
        
    end
    
    methods (Access = private, Static)
        
        function s = createScalarProductSettings()
            s.scalarProductSettings.femSettings = [];
            s.epsilon = [];
        end
        
    end
    
end