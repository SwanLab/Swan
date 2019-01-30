classdef testAmplificatorsPdependency < ...
        testShowingError & ...
        testLoadStoredVariable & ...
        testStoredComputedChecker
    
    
    properties (Access = protected)
        variablesToStore = {'P'};
        tol = 1e-6;
        testName = 'AmplificatorsPdependency.mat';
    end
    
    properties (Access = private)
        nExp
        pExp
        firstInvP
        m1 = 0.2;
        m2 = 0.8;
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
    end
    
    methods (Access = public)
        
        function obj = testAmplificatorsPdependency()
            obj.init();
            obj.createHomogenizers();
            obj.createShapeFunctions();
            obj.printSmoothAndNonSmoothVariables();            
            obj.computeFirstAmplicatorComponents();
            obj.selectComputedVar();
            obj.plotFirstAmplificatorComponent();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nExp = 6;
            obj.pExp = 2.^[1:obj.nExp];
        end
        
        function createHomogenizers(obj)
            obj.createSmoothRectInclusionHomogenizer();
            obj.createNonSmoothRectInclusionHomogenizer();
        end
        
        function createSmoothRectInclusionHomogenizer(obj)
            f    = strcat([obj.testName,obj.smoothStr]);
            p    = obj.printing;
            m1v  = obj.m1;
            m2v  = obj.m2;
            i    = obj.iter;
            obj.hSmooth = NumericalSmoothRectangleHomogenizer(f,p,m1v,m2v,i);
        end
        
        function createNonSmoothRectInclusionHomogenizer(obj)
            f    = strcat([obj.testName,obj.nonSmoothStr]);
            p    = obj.printing;
            m1v  = obj.m1;
            m2v  = obj.m2;
            i    = obj.iter;
            obj.hNonSmooth = NumericalRectangleHomogenizer(f,p,m1v,m2v,i);
        end
        
        function createShapeFunctions(obj)
            obj.shapeFunSmooth    = obj.createShapeFunction(obj.hSmooth);
            obj.shapeFunNonSmooth = obj.createShapeFunction(obj.hNonSmooth);
        end
        
        function sF = createShapeFunction(obj,homog)
            strain = obj.computeStrainWithCanonicalStress(homog);
            settings = homog.getSettings();
            sF = ShFunc_StressNorm(settings);
            sF.filter.preProcess();
            sF.setVstrain(strain);
        end
        
        function computeFirstAmplicatorComponents(obj)
            obj.firstInvP = ones(obj.nExp,2);
            for ip = 1:obj.nExp
                p = obj.pExp(ip);
                obj.firstInvP(ip,1) = obj.obtainFirstPinverseFromHomogenizer(obj.shapeFunSmooth,obj.hSmooth,p);
                obj.firstInvP(ip,2) = obj.obtainFirstPinverseFromHomogenizer(obj.shapeFunNonSmooth,obj.hNonSmooth,p);
            end
            
        end
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.firstInvP;
        end
        
        function plotFirstAmplificatorComponent(obj)
            figureID = figure(1);
            h{1} = plot(obj.pExp,obj.firstInvP,'-+');
            legend({obj.smoothStr,obj.nonSmoothStr})
            figureName = 'InverseFirstPcompwithPnorm';
            outPutFigName = ['/home/alex/Dropbox/Amplificators/Images/',figureName];
            xlabelName = 'Pnorm';
            fp = figurePlotter(figureID,h,xlabelName);
            fp.print(outPutFigName)
        end
        
        
        function f = obtainFirstPinverseFromHomogenizer(obj,sF,homog,p)
            sF.setPnorm(p);
            ls = homog.getLevelSet();
            sF.computeCostAndGradient(ls)
            Pv = sF.getValue();
            f = 1/((Pv)^(1/p));
        end
        
        function strain = computeStrainWithCanonicalStress(obj,homog)
            stress = [1,0,0]';
            Ch = homog.getCh();
            strain = Ch\stress;
        end
        
        function printSmoothAndNonSmoothVariables(obj)
            obj.printVariables(obj.shapeFunSmooth,obj.hSmooth,obj.smoothStr);
            obj.printVariables(obj.shapeFunNonSmooth,obj.hNonSmooth,obj.nonSmoothStr);
        end
        
        function printVariables(obj,sF,homog,strCase)
            p = 2;
            sF.setPnorm(p);
            ls = homog.getLevelSet();
            sF.computeCostAndGradient(ls)
            obj.createPostProcess(sF,strCase);
            obj.print(sF,homog);
        end
        
        function createPostProcess(obj,sF,strCase)
            obj.createPostProcessDataBase(sF,strCase);
            postCase = 'ElasticityMicroAndLevelSet';
            obj.postProcess = Postprocess(postCase,obj.dataBaseForPost);
        end
        
        function createPostProcessDataBase(obj,sh,strCase)
            microProb          = sh.getPhysicalProblems();
            dI.mesh            = microProb{1}.mesh;
            dI.outName         = strcat([obj.testName,strCase]);
            ps = PostProcessDataBaseCreatorWithNoGaussData(dI);
            obj.dataBaseForPost = ps.getValue();
        end
        
        function print(obj,sF,homog)
            d.dens     = homog.getDensity();
            d.levelSet = homog.getLevelSet();
            microProb  = sF.getPhysicalProblems();
            d.quad = microProb{1}.element.quadrature;
            d.microVar = microProb{1}.variables;
            obj.postProcess.print(obj.iter,d);
        end
        
    end
    
    
end