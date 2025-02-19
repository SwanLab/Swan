classdef ResultsPrinter < GiDPrinter
    
    properties (Access = protected)
        ext = 'res'
        fields
        headPrinter
    end
    
    properties (Access = protected, Abstract)
        simulationStr
        hasGaussData
    end
    
    methods (Access = public, Static)
        
        function rP = create(resultCase,d)
            f = ResultsPrinterFactory();
            rP = f.create(resultCase,d);
        end
        
    end
    
    methods (Access = public)
        
        function print(obj,iter,d)
            obj.iter = iter;
            obj.createFileName(iter);
            obj.openFile();
            obj.storeFieldsToPrint(d);
            dh = obj.createHeadPrinterCommonInfo();
            obj.createHeadPrinter(d,dh)
            obj.printHeader();
            obj.printResults(iter,obj.fileID);
            obj.closeFile();
        end
        
        function setSimulationStr(obj,s)
            obj.simulationStr = s;
        end
        
        function itHas = getHasGaussData(obj)
            itHas = obj.hasGaussData;
        end
        
        function h = getHeadPrinter(obj)
            h = obj.headPrinter;
        end
        
    end
    
    
    methods (Access = protected)
        
        function init(obj,d)
            obj.etype = d.etype;
            obj.ndim  = d.ndim;
            obj.outFileName = d.outFileName;
            obj.resultsDir = d.resultsDir;
        end
        
        function d = createScalarDataBase(obj,iter,fileID,fieldValues,fieldName, fieldPosition)
            d.fileID        = fileID;
            d.fieldValues   = fieldValues;
            d.fieldName     = fieldName;
            d.iter          = iter;
            d.fieldPosition = fieldPosition;
            d.simulationStr = obj.simulationStr;
        end
        
        function d = createVectorDataBase(obj,iter,fileID,fieldValues,fieldName,fieldPosition,fieldComponentName)
            d = obj.createScalarDataBase(iter,fileID,fieldValues,fieldName,fieldPosition);
            d.fieldComponentName = fieldComponentName;
        end
        
        function d = createScalarGaussDataBase(obj,varargin)
            d = obj.createScalarDataBase(varargin{:});
            d = obj.addGaussDescriptor(d);
        end
        
        function d = createVectorGaussDataBase(obj,varargin)
            d = obj.createVectorDataBase(varargin{:});
            d = obj.addGaussDescriptor(d);
        end
        
        function createHeadPrinter(obj,d,dh)
            if obj.getHasGaussData()
                obj.headPrinter = GaussHeadPrinter(d,dh);
            else
                obj.headPrinter = NoGaussHeadPrinter(dh);
            end
        end
        
    end
    
    methods (Access = private)
        
        function dh = createHeadPrinterCommonInfo(obj)
            dh.fileID = obj.fileID;
            dh.etype = obj.etype;
            dh.ndim = obj.ndim;
        end
        
        function d = addGaussDescriptor(obj,d)
            d.gaussDescriptor = 'Guass up?';
        end
        
        function printHeader(obj)
            obj.headPrinter.print();
        end
        
    end
    
    methods (Abstract, Access = protected)
        storeFieldsToPrint(obj)
    end
    
    methods (Abstract, Access = public)
        printResults(obj)
    end
    
end
