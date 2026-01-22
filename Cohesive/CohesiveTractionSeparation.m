classdef CohesiveTractionSeparation < handle
    
    properties (Access = public)
        lawType
        law
        normalCharLength
        tangencialCharLength
    end
    
    methods (Access = public)
        
        function obj = CohesiveTractionSeparation(cParams)
            obj.init(cParams) 
        end

        function t = evaluate(obj, disp)
            t = obj.law.evaluate(disp);
        end

        function d = derivative(obj,disp)
            d=obj.law.derivative(disp);
        end

        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.lawType = cParams.lawType;
            obj.normalCharLength = cParams.normalCharLength;
            obj.tangencialCharLength = cParams.tangencialCharLength;
            obj.sigmaMax = cParams.sigmaMax;
            
            s.normalCharLength = obj.normalCharLength;
            s.tangencialCharLength = obj.tangencialCharLength;
            s.sigmaMax = obj.sigmaMax;
            if isfield(cParams,'a'), s.a = cParams.a; end
            if isfield(cParams,'b'), s.b = cParams.b; end
        
            obj.law = CohesiveLawFactory.create(obj.lawType, s); 
              
        end
        

    end
    
end