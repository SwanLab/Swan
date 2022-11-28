classdef ContinuousDiscreteSectionsInterpolation < handle
    
    properties (Access = public)
       sectionArea
       sectionInertia
    end

    properties (Access = private)
        sectionAreaInformation
        sectionInertiaInformation
    end
       
    methods (Access = protected)

        function computeSectionAreaAndInertia(obj)
            obj.computeSectionArea(obj);
            obj.computeSectionInertia(obj);
        end
        
        function computeSectionArea(obj)
            s    = obj.designVariable;
            sNum = contToDiscrete(s);
            obj.sectionArea = interpolateSectionsInformation(sNum,obj.sectionAreaInformation);
        end

        function computeSectionInertia(obj)
            s    = obj.designVariable;
            sNum = contToDiscrete(s);
            obj.sectionInertia = interpolateSectionsInformation(sNum,obj.sectionInertiaInformation);
        end
        
        
    end
    
    methods (Access = protected, Static)

        function x = contToDiscrete(x)
            x = x*36 + 1;
            x = max(1,min(x,37));
        end

        function val = interpolateSectionsInformation(x,property)
            val    = zeros(length(s),1);
            intenger = fix(x);
            incr     = x - intenger;
            for i = 1:length(x)
                if intenger(i) == 37
                    val(i) = property(intenger(i));
                else
                    k1   = property(intenger(i));
                    k2   = property(intenger(i) + 1);
                    val(i) = (k1*(1-incr(i)) + k2*incr(i));
                end
            end
        end
        
    end
    
end