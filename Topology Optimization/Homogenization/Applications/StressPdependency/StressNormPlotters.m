classdef StressNormPlotters < handle
    
    properties (Access = private)
       plotters         
    end
    
    
    methods (Access = public)
       
        function obj = StressNormPlotters(d)
            obj.plotters{1} = StressAndMaxNormWithExponentPlotter(d);
            obj.plotters{2} = StressMaxWithMeshSizePlotter(d);                        
        end
        
        function plot(obj)
            obj.plotters{1}.plot();
            obj.plotters{2}.plot();
        end
                
    end
      
    
end




