classdef SettingsOptimizerUnconstrained < AbstractSettings
    
   properties (Access = protected)
      defaultParamsName = 'paramsOptimizerUnconstrained' 
   end
    
   properties (Access = public)
       nconstr
       target_parameters
       constraint_case
       
       epsilon
       scalarProductSettings
       lineSearchSettings
       
       e2
       filter
       printChangingFilter       
       filename
       ptype
   end
   
   methods (Access = public)

       function obj = SettingsOptimizerUnconstrained(varargin)
           if nargin == 1
            obj.loadParams(varargin{1});
           end
       end
       
   end
    
  
    
end