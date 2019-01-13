classdef ResultsPrinterFactory < handle
    
    properties (Access = private)
        dStandard
        printer
    end
    
    methods (Access = public)
        function p = create(obj,resultCase,d,dT)
            obj.createStandardDataBase(d.dStandard);
            obj.createPrinter(resultCase,d,dT);
            p = obj.printer;
        end
        
    end
    
    methods (Access = private)
        
        function createPrinter(obj,resultCase,d,dT)
            dS = obj.dStandard;           
            switch resultCase
                case 'Elasticity'
                    dG = d.dGauss;
                    p = ElasticityResultsPrinter(dS,dG);
                case 'ElasticityMicro'
                    dG = d.dGauss;
                    p = ElasticityMicroResultsPrinter(dS,dG);
                case 'TopOptProblem'                    
                   % dT = obj.createDataBaseForTopOpt(d.dStandard);
                    d.dStandard = dS;
                    p = TopOptResultsPrinter.create(d,dT,dS.hasGaussData);
                case 'DensityGauss'
                    dG = d.dGauss;
                    p = DensityGaussResultsPrinter(dS,dG);
                case 'LevelSet'
                    p = LevelSetResultsPrinter(dS);
                case 'Density'
                    p = DensityResultsPrinter(dS);
            end
            obj.printer = p;
        end
        
        function createStandardDataBase(obj,d)
            dS.etype           = d.etype;
            dS.ndim            = d.ndim;
            dS.ptype           = d.ptype;
            dS.outFileName     = d.outFileName;
            dS.hasGaussData    = d.hasGaussData;
            dS.resultsDir      = d.resultsDir;
            obj.dStandard = dS;
        end
        
    end
    
    methods (Access = private, Static)
        
        function dT = createDataBaseForTopOpt(d)
            dT.optimizer       = d.optimizer;
            dT.printMode       = d.printMode;
            dT.ShapeNames      = d.ShapeNames;
        end
    end
end
