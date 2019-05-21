classdef testPrintingPostProcessComputer < handle
    
    properties (Access = private)
        iter
        fileName
        postProcessCase
        postProcessCreateData
        postProcessPrintData
    end    
    
    methods (Access = public)
        
        function obj = testPrintingPostProcessComputer(d)
            obj.init(d);
            obj.loadPostProcessCreateData();
            obj.loadPostProcessPrintData();
            obj.print();
        end
        
        function print(obj)
            it = obj.iter;
            pC = obj.postProcessCase;
            cD = obj.postProcessCreateData;
            pD = obj.postProcessPrintData;
            postProcess = Postprocess(pC,cD);
            postProcess.print(it,pD);            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.iter            = d.iter;
            obj.fileName        = d.file;
            obj.postProcessCase = d.postCase;
        end
        
        function loadPostProcessCreateData(obj)
            loadData = load(obj.fileName);            
            cD = loadData.postProcessData.createData;
            obj.postProcessCreateData = cD;
        end
        
        function loadPostProcessPrintData(obj)
            loadData = load(obj.fileName);                        
            pD = loadData.postProcessData.printData;
            obj.postProcessPrintData = pD;
        end       
                
    end
    
end