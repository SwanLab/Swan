classdef ContinuousDiscreteSectionsInterpolation < MaterialInterpolation
    
    properties (Access = protected)
       dmu0
       dmu1
       dk0
       dk1
    end
       
    methods (Access = protected)

        function [A,I] = computeSectionAreaAndInertia(obj)
            A = computeSectionArea(obj);
            I = computeInertia(obj);
        end
        
        
        
    end
    
    methods (Access = protected, Static)
        
        function A = computeSectionArea(obj)
            s = obj.designVariable
        end
        
    end
    
end