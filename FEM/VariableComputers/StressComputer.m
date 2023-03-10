classdef StressComputer < handle
    
    properties (Access = private)
        C
        strain
    end
    
    methods (Access = public)
        
        function obj = StressComputer(cParams)
            obj.init(cParams)
        end

        function stress = compute(obj)
            stress = obj.computeStress();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.C      = cParams.C;
            obj.strain = cParams.strain;
        end

        function stress = computeStress(obj)
            Cmat  = obj.C;
            strn  = permute(obj.strain,[2 3 1]);
            strn2(:,1,:,:) = strn;
            stress=squeeze(pagemtimes(Cmat,strn2));
            stress = permute(stress, [3 1 2]);
        end
        
    end
    
end