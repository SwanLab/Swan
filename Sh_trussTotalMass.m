classdef Sh_trussTotalMass < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        interpolator
        barLength
        designVariable
    end
    
    methods (Access = public)
        
        function obj = Sh_trussTotalMass(cParams)
            obj.init(cParams);
        end
        
        
        function computeFunction(obj)
            interp = obj.interpolator;
            interp.computeSectionArea();
            obj.value = interp.sectionArea;
        end

        function computeGradient(obj)            
            xR = obj.designVariable(1:nVar);
            xT = obj.designVariable(nVar+1:end);
            gR = 4*pi*xT;
            gT = 4*pi*xR;
            obj.gradient = [gR;gT];         
        end
        
        
    end
    
    methods (Access = private)

        function obj = init(obj,cParams)
            obj.interpolator   = cParams.interp;
            obj.designVariable = cParams.designVar;
            obj.varN           = length(obj.designVariable)/2;
            obj.barLength      = load(cParams.barLength);
        end

    end
    
end