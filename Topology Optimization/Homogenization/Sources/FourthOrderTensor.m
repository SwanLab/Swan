classdef FourthOrderTensor < Tensor
    
    properties (Access = public)
        tensor
    end
    
    properties            
            tensorVoigt 
            tensorVoigtInPlaneStress
            InverseTensorVoigt
            InverseTensorVoigtInPlaneStress
    end
    
    properties (Access = private)
        IndexTransformer
    end
    
    methods
        
        function obj = FourthOrderTensor()
            obj.IndexTransformer = TensorVoigtIndexTransformer();
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

