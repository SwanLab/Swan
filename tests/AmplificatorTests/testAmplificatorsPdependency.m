classdef testAmplificatorsPdependency < ...
        testShowingError & ...
        testLoadStoredVariable & ...
        testStoredComputedChecker
    
    
    properties (Access = protected)
        variablesToStore = {'P'};
        tol = 1e-6;
        testName = 'AmplificatorsPdependency';
    end
    
    properties (Access = private)
        nExp
        pExp
        firstInvP
        m1 = 0.2;
        m2 = 0.8;
        pnorm = 4;
        printing = false;
        iter = 0;
        hSmooth
        hNonSmooth
        shapeFunSmooth
        shapeFunNonSmooth
        smoothStr = 'Smooth';
        nonSmoothStr = 'NonSmooth';
        postProcess
        stressShape
        dataBaseForPost
        outFileName
        nDB
        stressMaxSmooth
        stressMaxNonSmooth
        expName
    end
    
    methods (Access = public)
        
        function obj = testAmplificatorsPdependency()
            obj.init();
            experiments = {'001','0005','0001'};
            for iexp = 1:length(experiments)
                obj.expName = experiments{iexp};
                obj.createHomogenizers();
                obj.createShapeFunctions();
                obj.printSmoothAndNonSmoothVariables();
                obj.computeFirstAmplicatorComponents();
                obj.computeMaxStresses();
                obj.selectComputedVar();
                obj.plotFirstAmplificatorComponent();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nExp = 6;
            obj.pExp = 2.^(1:obj.nExp);
        end
        
        function createHomogenizers(obj)
            gSmooth    = ['/home/alex/Desktop/MeshElips/FreeFem/Alex',obj.expName,'.msh'];
            gNonSmooth = ['/home/alex/Desktop/MeshElips/FreeFem/AlexRect',obj.expName,'.msh'];
            obj.hSmooth    = obj.createHomogenizer(obj.smoothStr,gSmooth);
            obj.hNonSmooth = obj.createHomogenizer(obj.nonSmoothStr,gNonSmooth);
        end
        
        function homog = createHomogenizer(obj,strCase,gFile)
            microFile = strcat([obj.testName,strCase,obj.expName]);
            obj.createFemInputData(gFile,microFile);
            nHDB = NumericalHomogenizerDataBase([microFile,'.m']);
            obj.nDB = nHDB.dataBase;
            d = obj.nDB;
            d.outFileName = microFile;
            d.print = false;
            d.levelSetDataBase.levelSetType = 'full';
            homog = NumericalHomogenizer(d);
        end
        
        function createFemInputData(obj,gmsFile,outfileName)
            oD = fullfile('Output',outfileName);
            obj.outFileName = fullfile(oD,[outfileName,'.m']);
            c = GmsFile2FemMatOoFileConverter(gmsFile,oD,obj.outFileName);
            c.convert();
        end
        
        function createShapeFunctions(obj)
            obj.shapeFunSmooth    = obj.createShapeFunction(obj.hSmooth,obj.smoothStr);
            obj.shapeFunNonSmooth = obj.createShapeFunction(obj.hNonSmooth,obj.nonSmoothStr);
        end
        
        function sF = createShapeFunction(obj,homog,strCase)
            strain = obj.computeStrainWithCanonicalStress(homog);
            dB.filter = 'P1';
            dB.optimizer = 'SLERP';
            dB.pdim = '2D';
            microFile = strcat([obj.testName,strCase,obj.expName]);
            dB.filename = microFile;
            dB.ptype = 'MICRO';
            dB.TOL      = obj.nDB.materialDataBase.matProp;
            dB.material = obj.nDB.materialDataBase.materialType;
            dB.method   = obj.nDB.materialInterpDataBase.method;
            dB.pdim = obj.nDB.pdim;
            sF = ShFunc_StressNorm(dB);
            sF.filter.preProcess();
            sF.setVstrain(strain);
        end
        
        function computeFirstAmplicatorComponents(obj)
            obj.firstInvP = ones(obj.nExp,2);
            for ip = 1:obj.nExp
                p = obj.pExp(ip);
                obj.firstInvP(ip,1) = obj.obtainFirstPinverseFromHomogenizer(obj.shapeFunSmooth,p);
                obj.firstInvP(ip,2) = obj.obtainFirstPinverseFromHomogenizer(obj.shapeFunNonSmooth,p);
            end
            
        end
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.firstInvP;
        end
        
        function plotFirstAmplificatorComponent(obj)
            figureID = figure(1);
            h{1} = plot(obj.pExp,obj.firstInvP(:,1),'b-+',obj.pExp,obj.firstInvP(:,2),'r-+');
            hold on;
            sigMaxToPlot(:,1) = 1/obj.stressMaxSmooth*ones(length(obj.pExp),1);
            sigMaxToPlot(:,2) = 1/obj.stressMaxNonSmooth*ones(length(obj.pExp),1);
            h{2} = plot(obj.pExp,sigMaxToPlot(:,1),'b--',obj.pExp,sigMaxToPlot(:,2),'r--');
            legend({obj.smoothStr,obj.nonSmoothStr,'SmoothMaxNorm','NonSmoothMaxNorm'});
            figureName = 'InverseFirstPcompwithPnorm';
            outPutFigName = ['/home/alex/Dropbox/Amplificators/Images/',figureName];
            xlabelName = 'Pnorm';
            fp = figurePlotter(figureID,h,xlabelName);
            %fp.print(outPutFigName)
        end
        
        
        function f = obtainFirstPinverseFromHomogenizer(obj,sF,p)
            sF.setPnorm(p);
            Pv = sF.computeCostWithFullDomain();
            f = 1/((Pv)^(1/p));
        end
        
        function computeMaxStresses(obj)
            obj.stressMaxSmooth    = obj.computeMaxStress(obj.shapeFunSmooth);
            obj.stressMaxNonSmooth = obj.computeMaxStress(obj.shapeFunNonSmooth);
        end
        
        function s = computeMaxStress(obj,sF)
            s = sF.computeMaxStressWithFullDomain();
        end
        
        function printSmoothAndNonSmoothVariables(obj)
            obj.printVariables(obj.shapeFunSmooth,obj.smoothStr);
            obj.printVariables(obj.shapeFunNonSmooth,obj.nonSmoothStr);
        end
        
        function printVariables(obj,sF,strCase)
            p = 2;
            sF.setPnorm(p);
            sF.computeCostWithFullDomain();
            obj.createPostProcess(sF,strCase);
            obj.print(sF);
        end
        
        function createPostProcess(obj,sF,strCase)
            obj.createPostProcessDataBase(sF,strCase);
            postCase = 'ElasticityMicro';
            obj.postProcess = Postprocess(postCase,obj.dataBaseForPost);
        end
        
        function createPostProcessDataBase(obj,sh,strCase)
            microProb          = sh.getPhysicalProblems();
            dI.mesh            = microProb{1}.mesh;
            dI.outName         = strcat([obj.testName,strCase,obj.expName]);
            ps = PostProcessDataBaseCreator(dI);
            obj.dataBaseForPost = ps.getValue();
        end
        
        function print(obj,sF)
            microProb  = sF.getPhysicalProblems();
            d.quad     = microProb{1}.element.quadrature;
            d.fields   = microProb{1}.variables;
            obj.postProcess.print(obj.iter,d);
        end
        
    end
    
    
    methods (Access = private, Static)
        
        function strain = computeStrainWithCanonicalStress(homog)
            stress = [1,0,0]';
            Ch = homog.Ch();
            strain = Ch\stress;
        end
        
    end
    
    
end