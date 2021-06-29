classdef VoigtPlaneStressHomogHomogenizer < SequentialLaminateHomogenizer 
    
    properties (Access = private)

    end
    
    methods (Access = public)
        
        function obj = VoigtPlaneStressHomogHomogenizer(c0,c1,dir,mi,theta)
            obj.init(c0,c1,dir,mi,theta);
            obj.computeAnisotropicTensors();
            obj.transformTensorsToVoigt()
            obj.makeTensorsPlaneStress()
            obj.homogenize()
            obj.planeStressHomogenizedTensor = obj.homogTensor;
        end

    end
    
    methods (Access = private)

         function makeTensorsPlaneStress(obj)
            obj.makeFiberAndMatrixTensorsPlaneStress()
            obj.makeAnisotropicTensorsPlaneStress()
         end
         
         function makeFiberAndMatrixTensorsPlaneStress(obj)
             obj.matTen = PlaneStressTransformer.transform(obj.matTen);
             obj.fibTen = PlaneStressTransformer.transform(obj.fibTen);
         end
         
        function makeAnisotropicTensorsPlaneStress(obj)
            for irank = 1:obj.rank
                ani = obj.anisotropicTensors{irank};
                aniPS = PlaneStressTransformer.transform(ani);
                obj.anisotropicTensors{irank} = aniPS;
            end
        end
         
    end
    
    
end

