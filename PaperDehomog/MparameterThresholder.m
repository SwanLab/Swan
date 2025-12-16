classdef MparameterThresholder < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mV
    end
    
    properties (Access = private)
       minLengthInUnitCell
       m0V
       m1V
    end
    
    methods (Access = public)
        
        function obj = MparameterThresholder(cParams)
            obj.minLengthInUnitCell = cParams.minLengthInUnitCell;
            obj.m0V = 0;
            obj.m1V = 1;
        end
        
        function m = thresh(obj,m)  
            m0 = 0.001;min(m);%obj.m0V;
            m1 = 0.999;max(m);%obj.m1V;
            t = obj.minLengthInUnitCell;
            isAlmostClosed  = m <= t/2;
            isRoughlyClosed = (t/2 <= m) & (m <= t);           
            isRoughlyOpened = ((1-t) <= m) & (m <= (1-t/2));
            isAlmostOpened  = (1-t/2) <= m;
            m(isAlmostClosed)  = m0;
            m(isRoughlyClosed) = t(isRoughlyClosed);
            m(isRoughlyOpened) = 1-t(isRoughlyOpened);
            m(isAlmostOpened)  = m1;
        end
        
    end    
    
    methods (Access = private, Static)
        
        function itIs = isMcloseToMmin(m,mMin,mMax)
           itIs = (m - mMin) < (mMax - mMin)/2;
        end
        
    end
    
end