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
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.lawType = cParams.lawType;
            obj.normalCharLength = cParams.normalCharLength;
            obj.tangencialCharLength = cParams.tangencialCharLength;
            obj.sigmaMax = cParams.sigmaMax;
            
            switch obj.lawType
                case 'Cubic'
                    s.normalCharLength = obj.normalCharLength;
                    s.tangencialCharLength = obj.tangencialCharLength;
                    s.sigmaMax = obj.sigmaMax;
                    obj.law = CubicCohesiveLaw(s);
                case 'Bilinear'
                    s.a = cParams.a;
                    s.b = cParams.b;
                    obj.law = BilinearCohesiveLaw(s);
            end
              
        end
        

    end
    
end