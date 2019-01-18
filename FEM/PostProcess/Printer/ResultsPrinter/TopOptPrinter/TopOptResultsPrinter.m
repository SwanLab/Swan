classdef TopOptResultsPrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr = 'TopOpt';
        headPrinter
        headerData
        printers
        dResultsPrinter
    end
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = TopOptResultsPrinter(d)
            obj.init(d);
            f = TopOptResultsPrinterFactory(d);
            obj.printers = f.getPrinters();
            obj.outFileName = d.outFileName;
            obj.resultsDir = d.resultsDir;
            obj.setSimulationStrToPrinters();
            obj.createHeaderPrinter();
        end
        
        function storeResultsInfo(obj,d)
            
            if obj.hasGaussData()
                phyPr = d.cost.ShapeFuncs{1}.getPhysicalProblem();
                q = phyPr.element.quadrature;                
                obj.headerData.ngaus = q.ngaus;
                obj.headerData.posgp = q.posgp';
                obj.headerData.gaussDescriptor = 'Guass up?';
            end
            
            for iprinters = 1:numel(obj.printers)
                obj.printers{iprinters}.storeResultsInfo(d);
            end
        end
        
    end
    
    methods (Access = protected)
        
        function  printResults(obj)
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.print(obj.iter);
            end
        end
        
        function printHeader(obj)
            
            if obj.hasGaussData()
                obj.headerData.fileID = obj.fileID;
                obj.headerData.etype = obj.etype;
                obj.headerData.ndim  = obj.ndim;                
                obj.headPrinter.print(obj.headerData);
            else
                obj.headerData.fileID = obj.fileID;
                obj.headPrinter.print(obj.headerData);
            end
            

        end
        
    end
    
    methods (Access = private)
        
        function setSimulationStrToPrinters(obj)
            s = obj.simulationStr;
            for iprinters = 1:numel(obj.printers)
                obj.printers{iprinters}.setSimulationStr(s);
            end
        end
            
        function itHas = hasGaussData(obj)
            itHas = false;
            for iprinter = 1:numel(obj.printers)
                hasGaussDataPrinter = obj.printers{iprinter}.hasGaussData();
                itHas = hasGaussDataPrinter || itHas;
            end
        end
        
        function createHeaderPrinter(obj)
            if obj.hasGaussData()
                obj.headPrinter = GaussHeadPrinter;
            else
                obj.headPrinter = NoGaussHeadPrinter;
                obj.headerData  = [];
            end
        end
        
    end

    
end