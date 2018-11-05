classdef FourthOrderTensor < Tensor
    
    properties (Access = public)
            tensorVoigt 
            tensorVoigtInPlaneStress
            InverseTensorVoigt
            InverseTensorVoigtInPlaneStress
    end
        
    methods
        
        function obj = FourthOrderTensor()

        end
        
        function computeTensorVoigt(obj)
            obj.tensorVoigt = Tensor2VoigtConverter.convert(obj);
        end
        
        function computeTensorVoigtInPlaneStress(obj)
            obj.tensorVoigtInPlaneStress = PlaneStressTransformer.transform(obj.tensorVoigt);
        end
        
        
        function createRandomTensor(obj)
            obj.tensor = rand(3,3,3,3);  
            obj.MakeMajorAndMinorSymmetrization();
        end
        
        function MakeMajorAndMinorSymmetrization(obj)
            Symmetrizer = FourthOrderSymmetrizer;
            obj.tensor = Symmetrizer.symmetrize(obj.tensor);
        end
         

    end
    
    
end

