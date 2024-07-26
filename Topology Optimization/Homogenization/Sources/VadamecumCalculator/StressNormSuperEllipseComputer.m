classdef StressNormSuperEllipseComputer < handle
    
    properties (Access = private)
        microProblem
        fName
        outputFolder
    end
    
    properties (Access = private)
        mx
        my
        q
        phi
        pNorm
        print
        fileName
        hasToCaptureImage
        testName
        iter
        meshBackground
        mesh
        hMesh
    end
    
    methods (Access = public)
        
        function obj = StressNormSuperEllipseComputer(cParams)
            obj.init(cParams);
        end
        
        function sPnorm = compute(obj)
            obj.createMicroProblem()
            sPnorm = obj.computePstressNorm();
        end
        
        function var = computeCellVariables(obj)
            obj.compute();
            mProblem = obj.microProblem;
            var.Ctensor = mProblem.variables.varFromCh.Chomog;
            var.tstress = mProblem.variables.varFromCh.tstress;
            var.tstrain = mProblem.variables.varFromCh.tstrain;
            var.displ   = mProblem.variables.varFromCh.tdisp;
            var.volume  = mProblem.mesh.computeVolume();
            var.mesh    = mProblem.mesh;
            
            quad = Quadrature.set(var.mesh.geometryType);
            quad.computeQuadrature('CONSTANT');
            volume = var.mesh.computeDvolume(quad);
            var.integrationVar.dV    = volume;
            var.integrationVar.nstre = size(var.tstress,1);
            var.integrationVar.ngaus = size(var.tstress,2);
            var.integrationVar.geoVol = obj.meshBackground.computeVolume();
        end
        
        function printImage(obj)
            outputName = [obj.fileName,'Print'];
            if obj.print
                    s.mesh       = obj.mesh;
                    s.outPutName = outputName;
                    printer = SuperEllipsePrinter(s);
                    printer.print();
                if obj.hasToCaptureImage
                    printer.captureImage()
                end
            end
        end
        
        function printStress(obj)
            if obj.print
                microP = obj.microProblem;
                outputName = [obj.fileName,'Print'];
                dI.mesh    =  microP.mesh;
                dI.outName = outputName;
                dI.pdim    = '2D';
                dI.ptype   = 'MICRO';
                ps = PostProcessDataBaseCreator(dI);
                dB = ps.getValue();
                postCase = 'ElasticityMicro';
                postProcess = Postprocess(postCase,dB);
                d.fields = microP.variables;
                d.quad   = microP.element.quadrature;
                postProcess.print(obj.iter,d);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mx       = cParams.mx;
            obj.my       = cParams.my;
            obj.q        = cParams.q;
            obj.phi      = cParams.phi;
            obj.pNorm    = cParams.pNorm;
            obj.print    = cParams.print;
            obj.fileName = cParams.fileName;
            obj.iter     = cParams.iter;
            obj.hMesh    = cParams.hMesh;
            obj.hasToCaptureImage = cParams.hasToCaptureImage;
            obj.computeFileNameAndOutputFolder();
        end
        
        function computeFileNameAndOutputFolder(obj)
            obj.fName = [obj.fileName];
            obj.outputFolder = fullfile(pwd,'Output',obj.fName);
        end
        
        function sPnorm2 = computePstressNorm(obj)
            stress = [sin(obj.phi) cos(obj.phi) 0];
            p = obj.pNorm;
            sPnorm2 = obj.microProblem.computeStressPnorm(stress,p);
        end
        
        function ls = createLevelSet(obj)
            sM.coord  = obj.meshBackground.coord;
            sM.connec = obj.meshBackground.connec;
            s.mesh = Mesh_Total(sM);
            
            s.widthH = obj.mx;
            s.widthV = obj.my;
            s.pnorm  = obj.q;
            s.type = 'smoothRectangle';
            
            s.levelSetCreatorSettings = s;
            s.type = 'LevelSet';
            
            s.scalarProductSettings.epsilon = 1;
            levelSet = LevelSet(s);
            ls = levelSet.value;
        end
        
        function createBackgroundMesh(obj)
            run('RVE_Square_Triangle_FineFine')
            a.connec = connec(:, 2:end);
            a.coord  = coord(:, 2:3);
            m = Mesh.create(a);
            obj.testName = 'RVE_Square_Triangle_FineFine';
            %obj.testName = 'RVE_Square_Triangle_Fine';
            obj.meshBackground = m; 
        end
        
        function createMesh(obj)
            s.fileName = obj.fileName;
            s.levelSet = obj.createLevelSet();
            s.meshBackground = obj.meshBackground;
            s.hMesh = obj.hMesh;
            mCreator = MeshCreatorFromLevelSetWithMMG(s);
            obj.mesh = mCreator.create();
        end
        
        function createMicroProblem(obj)
            obj.createBackgroundMesh();
            obj.createMesh();
            femSolver = Elastic_Problem_Micro.create(obj.testName);
            femSolver.setMesh(obj.mesh);
            props.kappa = .75;
            props.mu    = .375;
            femSolver.setMatProps(props);
            obj.microProblem = femSolver;
        end
        
    end
    
    
end