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
        
        postData
        
        iMxIndex
        iMyIndex
        mxV
        myV
        
    end
    
    methods (Access = public)
        
        function obj = VademecumCalculator(d)
            obj.init(d);
            obj.computeVademecumData();
        end
        
        function d = getData(obj)
            d.postData     = obj.postData;
            d.postData.mxV = obj.mxV;
            d.postData.myV = obj.myV;
            d.outPutPath   = obj.outPutPath;
            d.fileName     = obj.fileName;
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
%             nMx = 5;
%             nMy = 5;
%             obj.mxV = linspace(0.2,0.8,nMx);
%             obj.myV = linspace(0.2,0.8,nMy);
        end
        
        function computeVademecumData(obj)
            nMx = length(obj.mxV);
            nMy = length(obj.myV);
            for imx = 1:nMx
                for imy = 1:nMy
                    obj.storeIndex(imx,imy);
                    obj.generateMeshFile();
                    obj.computeNumericalHomogenizer();
                    obj.obtainVolume();   
                    obj.obtainChi();
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
            obj.postData.Ctensor(:,:,imx,imy) = obj.homog.Ch;
        end
        
        function obtainPtensor(obj)
            P = obj.homog.Ptensor;
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;
            obj.postData.Ptensor(:,:,imx,imy) = P.getValue();
        end
        
        function obtainVolume(obj)
            v = obj.homog.volume;
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;
            obj.postData.volume(imx,imy) = v;
        end
        
        function obtainTxi(obj)
            v = obj.homog.volume;
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;
            obj.postData.volume(imx,imy) = v;
        end        
        
        function obtainPinvTensor(obj)
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;
            P = obj.postData.Ptensor(:,:,imx,imy);
            obj.postData.PinvTensor(:,:,imx,imy) = inv(P);
        end        
                
    end
    
end