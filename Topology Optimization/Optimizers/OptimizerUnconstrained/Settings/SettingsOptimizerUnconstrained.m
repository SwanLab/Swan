classdef SettingsOptimizerUnconstrained < AbstractSettings
    
   properties (Access = protected)
      defaultParamsName = 'paramsOptimizerUnconstrained.json' 
   end
    
   properties (Access = public)
       targetParameters
       designVariable
       lagrangian
       type
       
       convergenceVars
       
       epsilon
       scalarProductSettings
       lineSearchSettings
       
       e2    
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