classdef VademecumCalculator < handle
    
    properties (Access = private)
        linesRead
        
        fileName
        printingDir
        gmsFile
        outPutPath
        
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
            obj.printingDir      = fullfile(pwd,'Output',obj.fileName);                        
            obj.gmsFile          = [fullfile(obj.printingDir,obj.fileName),'.msh'];
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
            d.mxV             = obj.mxV(obj.iMxIndex);
            d.myV             = obj.myV(obj.iMyIndex);
            d.fileName        = obj.fileName;
            d.freeFemFileName = obj.fileName;
            d.printingDir = obj.printingDir;
            fG = FreeFemMeshGenerator(d);
            fG.generate();
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
        
        function obtainChi(obj)
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