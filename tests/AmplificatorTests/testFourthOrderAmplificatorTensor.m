classdef testFourthOrderAmplificatorTensor < testShowingError
    
    properties (Access = protected)
        tol = 5e-2;
        testName = 'testShapeStressWithAmplificator';
    end
    
    properties (Access = private)
        fileName   
        printingDir
        gmsFile        
    end
    
    
    methods (Access = public)
        
        function obj = testFourthOrderAmplificatorTensor()
            obj.init();
            obj.generateMeshFile();
            obj.createFemMatOoInputData();
            obj.createNumericalHomogenizer();
        end
        
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            obj.error = 0;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName    = 'FourthOrderAmplificator';
            obj.printingDir = fullfile(pwd,'Output',obj.fileName);  
            obj.gmsFile     = [fullfile(obj.printingDir,obj.fileName),'.msh'];            
        end
        
        function generateMeshFile(obj)
            d.mxV         = 0.8;
            d.myV         = 0.2;
            d.fileName    = obj.fileName;
            d.printingDir = obj.printingDir;
            d.freeFemFileName = 'SmoothRectangle';
            fG = FreeFemMeshGenerator(d);
            fG.generate();
        end       
        
        function createFemMatOoInputData(obj)
            oD = fullfile('Output',obj.fileName);
            fullOutFile = fullfile(oD,[obj.fileName,'.m']);
            c = GmsFile2FemMatOoFileConverter(obj.gmsFile,oD,fullOutFile);
            c.convert();
        end           
        
        function createNumericalHomogenizer(obj)
            defaultDB = NumericalHomogenizerDataBase([obj.fileName,'.m']);
            dB = defaultDB.dataBase;
            dB.outFileName                   = obj.fileName;
            dB.print                         = true;
            dB.levelSetDataBase.levelSetType = 'full';
            obj.homog = NumericalHomogenizer(dB);                                            
        end
         
    end
    
    
end