classdef NumericalHomogenizerCreatorFromGmsFile < handle
    
    properties (Access = private)
        gmsFile        
        outFile
        homog 
        homogDataBase     
        print
    end
   
    methods (Access = public)
        
        function obj = NumericalHomogenizerCreatorFromGmsFile(d)
            obj.init(d)
            obj.createFemMatOoInputData();
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
            obj.gmsFile = d.gmsFile;
            obj.outFile = d.outFile;
            obj.print   = d.print;
        end
                 
        function createFemMatOoInputData(obj)
            oD = fullfile('Output',obj.outFile);
            fullOutFile = fullfile(oD,[obj.outFile,'.m']);
            c = GmsFile2FemMatOoFileConverter(obj.gmsFile,oD,fullOutFile);
            c.convert();
        end        
        
        function createNumericalHomogenizerDataBase(obj)
            defaultDB = NumericalHomogenizerDataBase([obj.outFile,'.m']);
            dB = defaultDB.dataBase;
            dB.outFileName                   = obj.outFile;
            dB.print                         = obj.print;
            dB.levelSetDataBase.levelSetType = 'full';
            obj.homogDataBase = dB;
        end  
        
        function createNumericalHomogenizer(obj)
            d = obj.homogDataBase;
            obj.homog = NumericalHomogenizer(d);                                
        end
        
    end
    
    
end