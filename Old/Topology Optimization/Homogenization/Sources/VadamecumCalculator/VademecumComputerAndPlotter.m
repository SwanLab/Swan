classdef VademecumComputerAndPlotter < handle
    
    properties (GetAccess = public, SetAccess = private)
        vademecumData
    end
    
    properties (Access = private)
        fileName
        outPutPath
        vademecumPlotter
        vadVariables        
    end
    
    methods (Access = public)
        
        function obj = VademecumComputerAndPlotter(d)
            obj.init(d)
        end
        
        function compute(obj)
            obj.calculateVademecumPlotterData();
            obj.plotVademecum();
            obj.printData();
        end
        
        function fI = getFeasibleIndex(obj)
            vp = obj.vademecumPlotter;
            fI = vp.feasibleIndex;
        end
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.fileName   = d.fileName;
            obj.outPutPath = d.outPutPath;
        end
        
        function calculateVademecumPlotterData(obj)
            obj.calculateVademecumData();
            obj.transformVademecumCellDataInMatrix();
            d = obj.vademecumData;
            d.fileName   = obj.fileName;
            d.outPutPath = obj.outPutPath;
            obj.vademecumData = d;
        end
        
        function calculateVademecumData(obj)
            d.fileName = obj.fileName;
            d.outPutPath = obj.outPutPath;
            vc = VademecumPinvComputer(d);
            vc.compute();
            obj.vadVariables = vc.vademecumData;
        end
        
        function transformVademecumCellDataInMatrix(obj)
            mxV = obj.vadVariables.domVariables.mxV;
            myV = obj.vadVariables.domVariables.myV;            
            for imx = 1:length(mxV)
                for imy = 1:length(myV)
                    vad = obj.vadVariables.variables{imx,imy};
                    variables.C(:,:,imx,imy)     = vad.Ctensor;
                    variables.invP(:,:,imx,imy)  = vad.Pinvtensor;
                    variables.Ptensor(:,:,imx,imy) = vad.Ptensor;
                    variables.volume(imx,imy)    = vad.volume;
                end
            end            
            obj.vademecumData = variables;
            obj.vademecumData.mxV = mxV;
            obj.vademecumData.myV = myV;                        
        end        
        
        function plotVademecum(obj)
            d = obj.vademecumData;
            obj.vademecumPlotter = VademecumPlotters(d);
            obj.vademecumPlotter.plot();
        end
        
        function printData(obj)
            d.outPutPath = [obj.outPutPath,obj.fileName];
            d.Ctensor  = obj.vademecumData.C;
            d.Ptensor  = obj.vademecumData.Ptensor;
            d.volume   = obj.vademecumData.volume;
            p = VademecumDataPrinter(d);
           % p.print();
        end
        
    end
    
end