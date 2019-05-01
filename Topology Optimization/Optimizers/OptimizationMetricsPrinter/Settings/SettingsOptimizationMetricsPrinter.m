classdef SettingsOptimizationMetricsPrinter < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsOptimizationMetricsPrinter'
    end
    
    properties (Access = public)
      fileName
      shallPrint
      optimizer
      cost
      constraint
    end
    
     methods (Access = public)
        
        function obj = SettingsOptimizationMetricsPrinter(varargin)
            if nargin == 1
                    obj.loadParams(varargin{1});
            end
        end 
        
     end
     
end