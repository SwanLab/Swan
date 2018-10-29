classdef HomogVoigtPlaneStressHomogenizer < SequentialLaminateHomogenizer
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = HomogVoigtPlaneStressHomogenizer(c0,c1,dir,mi,theta)
            obj.init(c0,c1,dir,mi,theta);
            obj.computeAnisotropicTensors();
            obj.homogenize()            
            obj.transformToVoigt()
            obj.makeHomogenizedTensorPlaneStress();
        end               
    end
    
    methods (Access = protected)
        
         function transformToVoigt(obj)
            obj.transformHomogTensorInVoigt();
         end        
         
    end
    
    methods (Access = private)
        
        function transformHomogTensorInVoigt(obj)
            tens = FourthOrderTensor();
            tens.setValue(obj.homogTensor);
            obj.homogTensor = Tensor2VoigtConverter.convert(tens);            
        end
        
    end
    
end

