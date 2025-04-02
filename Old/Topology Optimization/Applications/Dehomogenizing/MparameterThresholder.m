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
            obj.init(cParams);
            obj.m0V = 0;
            obj.m1V = 1;
        end
        
        function m = thresh(obj,m)
            obj.mV = m;
            obj.thresholdingWhenTooSmallCell();
            obj.thresholdingLimitCases();
            m = obj.mV;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.minLengthInUnitCell = cParams.minLengthInUnitCell;
        end
        
        function itIs = isCellTooSmall(obj)
            t = obj.minLengthInUnitCell;
            itIs = t > 1/2;
        end
        
        function thresholdingWhenTooSmallCell(obj)
            isSmall = obj.isCellTooSmall();
            m = obj.mV;
            m0 = obj.m0V;
            m1 = obj.m1V;
            mT = m;
            mT(isSmall) = m0 + (m1-m0)*heaviside(m(isSmall)-1/2);
            obj.mV = mT;
        end
        
        function thresholdingLimitCases(obj)
            obj.thresholdingLeftCellFrame();
            obj.thresholdingRightCellFrame();
        end
        
        function thresholdingLeftCellFrame(obj)
           m = obj.mV; 
           t = obj.minLengthInUnitCell;
           mMin = obj.m0V*ones(size(t));
           mMax = mMin + t;  
           obj.mV = obj.thresholdCellFrame(m,mMin,mMax); 
        end
        
        function thresholdingRightCellFrame(obj)
           m = obj.mV;
           t = obj.minLengthInUnitCell;
           mMax = obj.m1V*ones(size(t));
           mMin = mMax - t;
           obj.mV = obj.thresholdCellFrame(m,mMin,mMax); 
        end
        
        function m = thresholdCellFrame(obj,m,mMin,mMax)
           isMin = obj.isMin(m,mMin,mMax);
           isMax = obj.isMax(m,mMin,mMax);
           m(isMin) = mMin(isMin);
           m(isMax) = mMax(isMax);
        end
        
        function itIs = isMin(obj,m,mMin,mMax)
           isFrame        = obj.isInFrameAndCellNotSmall(m,mMin,mMax);
           isCloseToMmin  = obj.isMcloseToMmin(m,mMin,mMax);
           itIs           = isFrame & isCloseToMmin;
        end
        
        function itIs = isMax(obj,m,mMin,mMax)
           isFrame        = obj.isInFrameAndCellNotSmall(m,mMin,mMax);
           isCloseToMmax  = ~obj.isMcloseToMmin(m,mMin,mMax);
           itIs           = isFrame & isCloseToMmax;
        end
        
        function itIs = isInFrameAndCellNotSmall(obj,m,mMin,mMax)
           isInFrame      = mMin < m & m < mMax;
           isCellNotSmall = ~obj.isCellTooSmall();
           itIs = isInFrame & isCellNotSmall;
        end
        
    end
    
    methods (Access = private, Static)
        
        function itIs = isMcloseToMmin(m,mMin,mMax)
           itIs = (m - mMin) < (mMax - mMin)/2;
        end
        
    end
    
end