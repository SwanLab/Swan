classdef VademecumCalculator < handle
    
    properties (Access = private)
        linesRead
        
        fileName
        printingDir
        filePath
        gmsFile
        freeFemFile
        freeFemModelFile
        outPutPath
        
        inputData
        homog
        Ctensor
        Ptensor
        PinvTensor
        
        iMxIndex
        iMyIndex
        mxV
        myV
    end
    
    methods (Access = public)
        
        function obj = VademecumCalculator(d)
            obj.init(d);
            obj.computeVademecumData();
            obj.makePlots();
            obj.printData();
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj,d)
            obj.fileName         = d.fileName;
            obj.outPutPath       = d.outPutPath;
            obj.freeFemModelFile = fullfile('Input',obj.fileName,[obj.fileName,'Model','.edp']);
            obj.printingDir      = fullfile(pwd,'Output',obj.fileName);
            obj.filePath         = fullfile(obj.printingDir,obj.fileName);
            obj.gmsFile          = [obj.filePath,'.msh'];
            obj.freeFemFile      = [obj.filePath,'.edp'];
            nMx = 20;
            nMy = 20;
            obj.mxV = linspace(0.01,0.99,nMx);
            obj.myV = linspace(0.01,0.99,nMy);
        end
        
        function computeVademecumData(obj)  
            nMx = length(obj.mxV);
            nMy = length(obj.myV);
            for imx = 1:nMx
                for imy = 1:nMy
                    obj.storeIndex(imx,imy);
                    obj.generateMeshFile();
                    obj.computeNumericalHomogenizer();
                    obj.obtainTensors();
                    (imy + nMx*(imx -1))/(nMx*nMy)*100
                end
            end                        
        end
        
        function storeIndex(obj,imx,imy)
            obj.iMxIndex = imx;
            obj.iMyIndex = imy;
        end
        
        function generateMeshFile(obj)
            obj.createInputData();
            obj.printFreeFemFile();
            obj.computeMeshWithFreeFem();
        end
        
        function createInputData(obj)
            obj.inputData.reals   = obj.createRealInputData();
            obj.inputData.strings = obj.createStringInputData();
        end
        
        function dB = createRealInputData(obj)
            mx = obj.mxV(obj.iMxIndex);
            my = obj.myV(obj.iMyIndex);
            dB{1,:} = {'p',4};
            dB{2,:} = {'mx',mx/2};
            dB{3,:} = {'my',my/2};
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
        
        function computeNumericalHomogenizer(obj)
            d.gmsFile = obj.gmsFile;
            d.outFile = obj.fileName;
            d.print   = false;
            nH = NumericalHomogenizerCreatorFromGmsFile(d);
            obj.homog = nH.getHomogenizer();
        end
        
        function obtainTensors(obj)
            obj.obtainCtensor();
            obj.obtainPtensor();
            obj.obtainPinvTensor();
        end
        
        function obtainCtensor(obj)
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;
            obj.Ctensor(:,:,imx,imy) = obj.homog.Ch;
        end
        
        function obtainPtensor(obj)
            P = obj.homog.Ptensor;
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;
            obj.Ptensor(:,:,imx,imy) = P.getValue();
        end
        
        function obtainPinvTensor(obj)
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;
            P = obj.Ptensor(:,:,imx,imy);
            obj.PinvTensor(:,:,imx,imy) = inv(P);
        end
        
        function makePlots(obj)
            d.mxV        = obj.mxV;
            d.myV        = obj.myV;
            d.C          = obj.Ctensor;
            d.invP       = obj.PinvTensor;
            d.hasToPrint = true; 
            d.outPutPath = obj.outPutPath;
            d.microName  = obj.fileName;
            p = VademecumPlotter(d);
            p.plot();
        end  
        
        function printData(obj)
            d.outPutPath = [obj.outPutPath,obj.fileName];            
            d.Ctensor  = obj.Ctensor;
            d.Ptensor  = obj.Ptensor;
            p = VademecumDataPrinter(d);
            p.print();            
        end
                   
        
    end
    
end