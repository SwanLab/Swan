classdef SettingsOptimizerUnconstrained < AbstractSettings
    
   properties (Access = protected)
      defaultParamsName = 'paramsOptimizerUnconstrained' 
   end
    
   properties (Access = public)
       target_parameters
       
       epsilon
       scalarProductSettings
       lineSearchSettings
       
       e2
       filter
       printChangingFilter       
       filename
       ptype
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