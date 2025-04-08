classdef VademecumVariablesComputer < handle
    
    properties (Access = public)
        vademecumData
    end
    
    properties (Access = protected)
        amplificatorInput
    end
    
    properties (Access = private)
        fileName
        samplingData
    end
    
    methods (Access = public)
        
        function compute(obj)
            obj.computeCellVariables();
            obj.computeAmplificatorComponents();
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,d)
            obj.fileName   = d.fileName;
        end
        
        function createAmplificatorInput(obj)
            obj.obtainIntegrationVariables();
            obj.obtainCellVariables();
        end
        
    end
    
    methods (Access = private)
        
        function computeCellVariables(obj)
            a  = load(obj.fileName);
            obj.vademecumData = a.d;
        end
        
        function computeAmplificatorComponents(obj)
            mxV = obj.vademecumData.domVariables.mxV;
            myV = obj.vademecumData.domVariables.myV;
            for imx = 1:length(mxV)
                for imy = 1:length(myV)
                    obj.specifySamplingCase(imx,imy);
                    obj.computeAmplificators(imx,imy);
                end
            end
        end
        
        function specifySamplingCase(obj,imx,imy)
            sD = obj.vademecumData.variables{imx,imy};
            obj.samplingData = sD;
        end
        
        function d = obtainIntegrationVariables(obj)
            intData   = obj.samplingData.integrationVar;
            d.nstre   = intData.nstre;
            d.V       = intData.geoVol;
            d.ngaus   = intData.ngaus;
            d.dV      = intData.dV;
            obj.amplificatorInput = d;
        end
        
        function d = obtainCellVariables(obj)
            d = obj.amplificatorInput;
            d.Ch      = obj.samplingData.Ctensor;
            d.tstress = obj.samplingData.tstress;
            obj.amplificatorInput = d;
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeAmplificators(obj)
    end
    
end