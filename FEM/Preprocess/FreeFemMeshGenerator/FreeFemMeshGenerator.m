classdef FreeFemMeshGenerator < handle
    
    properties (Access = private)
        mxV
        myV
        hMax
        intElements
        borderElements
        qNorm
        fileName
        printingDir
        filePath
        freeFemModelFile
        freeFemFileName
        freeFemFile
        inputData
        linesRead
        loadedFile
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
            obj.loadSettignsParams(d);
            obj.filePath         = fullfile(obj.printingDir,obj.fileName);
            obj.freeFemModelFile = fullfile('Input',obj.freeFemFileName,[obj.freeFemFileName,'Model','.edp']);
            obj.freeFemFile      = [obj.filePath,'.edp'];
        end
        
        function loadSettignsParams(obj,d)
          fields = fieldnames(d);
          for i = 1:length(fields)
                param = fields(i);
                param = param{1};
                if isprop(obj,param)
                    obj.(param) = d.(param);
                else
                    obj.warnOfInvalidCustomParams(param);
                end
          end
        end
        
        function createInputData(obj)
            obj.inputData.reals   = obj.createRealInputData();
            obj.inputData.strings = obj.createStringInputData();
        end
        
        function dB = createRealInputData(obj)
            mx = obj.mxV;
            my = obj.myV;
            dB{1,:} = {'p',obj.qNorm};
            dB{2,:} = {'mx',mx};
            dB{3,:} = {'my',my};
            dB{4,:} = {'Hmax',obj.hMax};
            dB{5,:} = {'elByInt',obj.intElements};
            dB{6,:} = {'elByBor',obj.borderElements};
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
        
         function warnOfInvalidCustomParams(obj,param)
            warning([param ' is not a property of ' class(obj)]);
        end
        
    end
    
end