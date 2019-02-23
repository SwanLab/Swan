classdef FreeFemMeshGenerator < handle
    
    properties (Access = private)
        mxV
        myV        
        fileName        
        printingDir
        filePath        
        freeFemModelFile
        freeFemFileName
        freeFemFile
        inputData
        linesRead        
    end
    
    methods (Access = public)
        
        function obj = FreeFemMeshGenerator(d)
           obj.init(d)             
        end
        
        function generate(obj)
            obj.createInputData();
            obj.printFreeFemFile();
            obj.computeMeshWithFreeFem();            
        end        
        
    end
        
    methods (Access = private)
        
        function init(obj,d)
            obj.mxV = d.mxV;
            obj.myV = d.myV;
            obj.fileName         = d.fileName;
            obj.printingDir      = d.printingDir;
            obj.freeFemFileName  = d.freeFemFileName;
            obj.filePath         = fullfile(obj.printingDir,obj.fileName);
            obj.freeFemModelFile = fullfile('Input',obj.freeFemFileName,[obj.freeFemFileName,'Model','.edp']);
            obj.freeFemFile      = [obj.filePath,'.edp'];                        
        end
               
        function createInputData(obj)
            obj.inputData.reals   = obj.createRealInputData();
            obj.inputData.strings = obj.createStringInputData();
        end
        
        function dB = createRealInputData(obj)
            mx = obj.mxV;
            my = obj.myV;
            dB{1,:} = {'p',4};
            dB{2,:} = {'mx',mx};
            dB{3,:} = {'my',my};
            dB{4,:} = {'Hmax',0.02};
            dB{5,:} = {'elByInt',20};
            dB{6,:} = {'elByBor',20};
        end
        
        function dB = createStringInputData(obj)
            dB{1,:} = {'OutName',obj.filePath};
        end
        
        function printFreeFemFile(obj)
            obj.readFile();
            obj.printFile();
        end
        
        function readFile(obj)
            fR = FreeFemFileReader(obj.freeFemModelFile);
            fR.read();
            obj.linesRead = fR.getDataBase();
        end
        
        function printFile(obj)
            d.fileName     = obj.fileName;
            d.printingDir  = obj.printingDir;
            d.linesToPrint = obj.linesRead;
            d.type = 'InputChange';
            fp = FreeFemFilePrinter.create(d);
            fp.setInputData(obj.inputData);
            fp.print();
        end
        
        function computeMeshWithFreeFem(obj)
            [~,~] = system(['FreeFem++ ',obj.freeFemFile]);
        end
        
    end
    
end