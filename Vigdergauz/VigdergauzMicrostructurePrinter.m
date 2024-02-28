classdef VigdergauzMicrostructurePrinter < handle
    
    properties (Access = private)
        settings
        unfittedMesh
        mesh
        jointMesh
        outputName
        inputFile
        pdim
        postProcess
        volumeMicro
        superEllipseRatio
        volumeMicroV
        superEllipseRatioV
        numericalHomogenizerDataBase
    end
    
    methods (Access = public)
        
        function obj = VigdergauzMicrostructurePrinter()
            obj.init();
            obj.createSamples();
            obj.printSamples();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.inputFile = 'VigergauzLevelSetInput';
            obj.outputName = 'VigdergauzPrinting';
            obj.pdim = '2D';
        end
        
        function createSamples(obj)
            obj.createVolumeSamples();
            obj.createSuperEllipseSamples();
        end
        
        function createVolumeSamples(obj)
            obj.volumeMicroV = [0.6,0.6,0.6,0.6,...
                                0.9 0.8 0.6,...
                                0.4 0.2 0.02];
        end
        
        function createSuperEllipseSamples(obj)
            obj.superEllipseRatioV = tan([pi/2-pi/8,pi/2-pi/7,pi/2-pi/6,pi/2-pi/5,...
                                      pi/4,pi/4,pi/4,...
                                      pi/4,pi/4,pi/4]);
            obj.checkFeasibility();
        end
        
        function printSamples(obj)
            for isample = 1:length(obj.superEllipseRatioV)
                obj.superEllipseRatio = obj.superEllipseRatioV(isample);
                obj.volumeMicro       = obj.volumeMicroV(isample);
                obj.createSettings();
                obj.createMesh();
                obj.createUnfittedMesh();
                obj.createJointMesh();
                obj.createPostProcess();
                obj.print(isample);                
            end
            
        end
        
        function createSettings(obj)
            setting = Settings(obj.inputFile);
            translator = SettingsTranslator();
            translator.translate(setting);
            fileName = translator.fileName;
            settingsTopOpt = SettingsTopOptProblem(fileName);
            obj.settings = settingsTopOpt;
            lsS = obj.settings.designVarSettings.levelSetCreatorSettings;
            lsS.vigdergauzDataBase.superEllipseRatio = obj.superEllipseRatio;
            lsS.vigdergauzDataBase.volumeMicro       = obj.volumeMicro;
            obj.settings.designVarSettings.levelSetCreatorSettings = lsS;
        end
        
        function createMesh(obj)
            s.coord = obj.settings.designVarSettings.femData.coord;
            s.connec = obj.settings.designVarSettings.femData.connec;
            obj.mesh = Mesh_Total(s);
        end
        
        function createUnfittedMesh(obj)
            s = obj.settings.designVarSettings;
            s.mesh = obj.mesh;
            s.scalarProductSettings.epsilon = 1e-3;
            designVariable = DesignVariable.create(s);
            obj.unfittedMesh = designVariable.getUnfittedMesh();
        end
        
        function createJointMesh(obj)
            coordInner     = obj.unfittedMesh.innerMesh.coord;
            connecInner    = obj.unfittedMesh.innerMesh.connec;
            coordCutInner  = obj.unfittedMesh.innerCutMesh.coord;
            connecCutInner = obj.unfittedMesh.innerCutMesh.connec;
            ncoord = size(coordInner,1);
            connecCutInner = connecCutInner + ncoord;
            coord = [coordInner;coordCutInner];
            connec = [connecInner;connecCutInner];
            obj.jointMesh = Mesh.create();
            obj.jointMesh.create(coord,connec);
        end
        
        function createPostProcess(obj)
            outName = obj.outputName;
            dI.mesh            = obj.jointMesh;
            dI.outName         = outName;
            dI.pdim            = '2D';
            dI.ptype           = 'MICRO';
            ps = PostProcessDataBaseCreator(dI);
            dB = ps.getValue();
            dB.printers = 'Density';
            postCase = 'Density';
            obj.postProcess = Postprocess(postCase,dB);
        end
        
        function print(obj,isample)
            d.x = ones(size(obj.jointMesh.coord(:,1),1),1);
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
        
        function checkFeasibility(obj)
            for isamples = 1:length(obj.volumeMicroV)
                if obj.isNotFeasible(isamples)
                    error('IsNotFeasible')
                end
            end
        end
        
        function itIsNot = isNotFeasible(obj,isample)
            isFirst  = obj.isFirstCondSatisfied(isample);
            isSecond = obj.isFirstCondSatisfied(isample);
            itIs = isFirst && isSecond;
            itIsNot = ~itIs;
        end
        
        function itIs = isFirstCondSatisfied(obj,isample)
            mxMax = 0.99;
            qMax = 10^6;
            rho = obj.volumeMicroV(isample);
            txi = atan(obj.superEllipseRatioV(isample));
            sE = SuperEllipseParamsRelator();            
            txiMax = sE.txiFromMxAndRho(mxMax,rho,qMax);            
            itIs = txi <= txiMax;
        end
        
        function itIs = isSecondCondSatisfied(obj,isample)
            myMax = 0.99;
            qMax = 10^6;
            rho = obj.volumeMicroV(isample);
            txi = atan(obj.superEllipseRatioV(isample));
            sE = SuperEllipseParamsRelator();            
            txiMin = sE.txiFromMyAndRho(myMax,rho,qMax);
            itIs = txiMin <= txi;
        end
        
    end
    
end
