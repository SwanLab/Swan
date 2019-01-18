classdef TopOptComplianceAndAdjointPrinter < TopOptShapePrinter
    
    methods (Access = public)
        
        function obj = TopOptComplianceAndAdjointPrinter(d)
            obj.createPrinters(d)
        end
        
        function itHas = hasGaussData(obj)
            itHas = true;
        end
        
        function storeResultsInfo(obj,shape)
            obj.storeComplianceResultInfo(shape);
            obj.storeAdjointResultInfo(shape);                        
        end
        
    end
    
    methods (Access = private)
        
        function storeComplianceResultInfo(obj,shape)
            phyPr = shape.getPhysicalProblem();
            d = obj.obtainVariablesAndQuad(phyPr);
            obj.printers{1}.storeResultsInfo(d);
        end
        
        function storeAdjointResultInfo(obj,shape)
            adjPr = shape.getAdjointProblem();
            d = obj.obtainVariablesAndQuad(adjPr);
            obj.printers{2}.storeResultsInfo(d);
        end
                
        function createPrinters(obj,d)
            obj.printers{1} = obj.createCompliancePrinter(d);
            obj.printers{2} = obj.createAdjointPrinter(d);
        end
        
        function p = createCompliancePrinter(obj,d)
            p =  ResultsPrinter.create('Elasticity',d);
            p.setStrVariablesNames('Stress','Strain','Disp');
        end
        
        function p = createAdjointPrinter(obj,d)
            p =  ResultsPrinter.create('Elasticity',d);
            p.setStrVariablesNames('StressAdj','StrainAdj','DispAdj');
        end
    end
end