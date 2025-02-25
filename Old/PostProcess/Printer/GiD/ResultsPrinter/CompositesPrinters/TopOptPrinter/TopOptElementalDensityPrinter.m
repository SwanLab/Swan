classdef TopOptElementalDensityPrinter < CompositeResultsPrinter
    
    methods (Access = public)
        
        function obj = TopOptElementalDensityPrinter(d)
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
        
        function createPrinters(obj,d)
            p = ResultsPrinter.create('DensityGauss',d);
            obj.printers{1} = p;
        end
        
        function storeFieldsToPrint(obj,d)
            phyPr = d.cost.shapeFunctions{1}.getPhysicalProblems();
            d.fields = obj.computeElementalDensity(d.x,phyPr{1});
            obj.printers{1}.storeFieldsToPrint(d);
        end
        
        function createHeadPrinter(obj,d,dh)
            phyPr = d.cost.shapeFunctions{1}.getPhysicalProblems();
            d.quad = phyPr{1}.element.quadrature;
            obj.printers{1}.createHeadPrinter(d,dh);
            h = obj.printers{1}.getHeadPrinter();
            obj.headPrinter = h;
        end
    end
    
    
    methods (Access = private, Static)
        
        function d = computeElementalDensity(x,phyPr)
            type = 'ElementalDensityCreatorByLevelSet'; %Only examples with ls by the moment
            d.levelSet = x;
            d.filterDataBase.shape = phyPr.element.interpolation_u.shape;
            d.filterDataBase.conec = phyPr.geometry.interpolation.T;
            d.filterDataBase.quadr = phyPr.element.quadrature;
            edc = ElementalDensityCreator.create(type,d);
            d = edc.getDensity();
        end
        
    end
    
end