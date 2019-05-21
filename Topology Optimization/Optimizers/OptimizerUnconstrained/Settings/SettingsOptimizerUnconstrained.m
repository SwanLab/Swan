classdef SettingsOptimizerUnconstrained < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsOptimizerUnconstrained.json'
    end
    
    properties (Access = public)
        type
        targetParameters
        designVariable
        lagrangian
        
        convergenceVars
        
        epsilon
        scalarProductSettings
        lineSearchSettings
        
        e2
        ub
        lb
    end
    
    methods (Access = public)
        
        function obj = SettingsOptimizerUnconstrained(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
end