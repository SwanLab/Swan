classdef VoigtHomogPlaneStressHomogenizer < SequentialLaminateHomogenizer 
    
    properties (Access = private)

    end
    
    methods (Access = public)
        
        function obj = VoigtHomogPlaneStressHomogenizer(c0,c1,dir,mi,theta)
            obj.init(c0,c1,dir,mi,theta);
            obj.computeAnisotropicTensors();
            obj.transformTensorsToVoigt()
            obj.homogenize()
            obj.makeHomogenizedTensorPlaneStress();
        end

    end
    
  
    
    
end

