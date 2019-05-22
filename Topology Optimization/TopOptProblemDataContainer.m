classdef TopOptProblemDataContainer < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = ''
    end
    
    properties (Access = public)
        caseFileName
        femFileName
        scale
        pdim
        nelem
        costFunctions
        costWeights
        constraintFunctions
        nConstraints
    end
    
    methods (Access = public)
        
        function obj = TopOptProblemDataContainer(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
end