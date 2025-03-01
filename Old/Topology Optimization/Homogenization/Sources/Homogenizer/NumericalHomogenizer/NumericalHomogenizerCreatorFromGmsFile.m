classdef NumericalHomogenizerCreatorFromGmsFile < handle
    
    properties (Access = private)
        gmsFile        
        outFile
        homog 
        homogDataBase     
        print
        iter
        hasToCaptureImage
    end
   
    methods (Access = public)
        
        function obj = NumericalHomogenizerCreatorFromGmsFile(d)
            obj.init(d)
            obj.createSwanInputData();
            obj.createNumericalHomogenizerDataBase();
            obj.createNumericalHomogenizer();
        end
        
        function h = getHomogenizer(obj)
            h = obj.homog;
        end
        
        function h = getHomogenizerDataBase(obj)
            h = obj.homogDataBase;
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.gmsFile           = d.gmsFile;
            obj.outFile           = d.outFile;
            obj.print             = d.print;
            obj.iter              = d.iter;
            obj.hasToCaptureImage = d.hasToCaptureImage;
        end
                 
        function createSwanInputData(obj)
            oD = fullfile('Output',obj.outFile);
            fullOutFile = fullfile(oD,[obj.outFile,'.m']);
            c = GmsFile2SwanFileConverter(obj.gmsFile,oD,fullOutFile);
            c.convert();
        end        
        
        function createNumericalHomogenizerDataBase(obj)
            defaultDB = NumericalHomogenizerDataBase([obj.outFile,'.m']);
            dB = defaultDB.dataBase;
            dB.outFileName                   = obj.outFile;
            dB.print                         = obj.print;
            dB.levelSetDataBase.type = 'full';
            dB.hasToCaptureImage = obj.hasToCaptureImage;
            obj.homogDataBase = dB;
        end  
        
        function createNumericalHomogenizer(obj)
            d = obj.homogDataBase;
            obj.homog = NumericalHomogenizer(d);
            obj.homog.iter = obj.iter;
            obj.homog.compute();                                
        end
        
    end
    
    
end