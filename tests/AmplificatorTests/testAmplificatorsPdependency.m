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
    end
    
    methods (Access = public)
        
        function obj = testAmplificatorsPdependency() 
            obj.init();
            obj.createSmoothRectInclusionHomogenizer();
            obj.createNonSmoothRectInclusionHomogenizer();
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
        
        function computeFirstAmplicatorComponents(obj)
            obj.firstInvP = ones(obj.nExp,2);
            for ip = 1:obj.nExp
               p = obj.pExp(ip);
               obj.firstInvP(ip,1) = obj.obtainFirstPinverseFromHomogenizer(obj.hSmooth,p);
               obj.firstInvP(ip,2) = obj.obtainFirstPinverseFromHomogenizer(obj.hNonSmooth,p);
            end

        end
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.firstInvP;
        end
        
        function plotFirstAmplificatorComponent(obj)
            figureID = figure(1);
            h{1} = plot(obj.pExp,obj.firstInvP,'-+');
            legend({'Smooth','NonSmooth'})
            figureName = 'InverseFirstPcompwithPnorm';
            outPutFigName = ['/home/alex/Dropbox/Amplificators/Images/',figureName];
            xlabelName = 'Pnorm';
            fp = figurePlotter(figureID,h,xlabelName);
            fp.print(outPutFigName)
        end
        
        function createSmoothRectInclusionHomogenizer(obj)
            f    = strcat(obj.testName);
            p    = obj.printing;
            m1v  = obj.m1;
            m2v  = obj.m2; 
            i    = obj.iter;
            obj.hSmooth = NumericalSmoothRectangleHomogenizer(f,p,m1v,m2v,i);                
        end
        
        function createNonSmoothRectInclusionHomogenizer(obj)
            f    = strcat(obj.testName);
            p    = obj.printing;
            m1v  = obj.m1;
            m2v  = obj.m2; 
            i    = obj.iter;
            obj.hNonSmooth = NumericalRectangleHomogenizer(f,p,m1v,m2v,i);                
        end
        
        function f = obtainFirstPinverseFromHomogenizer(obj,homog,p)
            strain = obj.computeStrainWithCanonicalStress(homog);
            settings = homog.getSettings();
            ls = homog.getLevelSet();
            sF = ShFunc_StressNorm(settings);
            sF.filter.preProcess();
            sF.setVstrain(strain);
            sF.setPnorm(p);            
            sF.computeCostAndGradient(ls)
            Pv = sF.getValue();
            f = 1/((Pv)^(1/p));
        end        
        
        function strain = computeStrainWithCanonicalStress(obj,homog)
            stress = [1,0,0]';
            Ch = homog.getCh();           
            strain = Ch\stress;
        end
        
    end
    
    
end