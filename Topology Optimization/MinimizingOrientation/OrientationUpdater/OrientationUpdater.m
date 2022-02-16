classdef OrientationUpdater < handle
    
    properties (GetAccess = public, SetAccess = protected)
        alpha
    end
    
    properties (Access = protected)
        eigenVectors
        eigenValues
        optimalIndexOrientation
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
           f = OrientationUpdaterFactory();
           obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
           
        function compute(obj,cParams)
            obj.eigenVectors = cParams.pD;
            obj.eigenValues  = cParams.pS;
            obj.computeOptimalIndexOrientation();
            obj.computeOrientation();
        end
        
    end
        
    methods (Access = private)
        
        function computeOrientation(obj)
            pD  = obj.eigenVectors;
            ind = obj.optimalIndexOrientation;
            isFirstOptimal  = ind == 1;
            isSecondOptimal = ind == 2;
            dir(:,isFirstOptimal)  = squeeze(pD(:,1,isFirstOptimal));
            dir(:,isSecondOptimal) = squeeze(pD(:,2,isSecondOptimal));
            obj.alpha = dir;
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeOptimalIndexOrientation(obj)
    end
    
end