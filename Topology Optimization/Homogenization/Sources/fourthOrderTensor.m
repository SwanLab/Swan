classdef fourthOrderTensor < handle
    
    
    properties
            tensor
            tensorVoigt 
            tensorVoigtInPlaneStress
            InverseTensorVoigt
            InverseTensorVoigtInPlaneStress
    end
    
    properties (Access = private)
        IndexTransformer
    end
    
    methods
        
        function obj = fourthOrderTensor()
            obj.IndexTransformer = TensorVoigtIndexTransformer();
        end
        
        function computeTensorVoigt(obj)
            obj.tensorVoigt = obj.RepresentTensorInVoigt(obj.tensor);
        end
        
        function computeTensorVoigtInPlaneStress(obj)
            obj.tensorVoigtInPlaneStress = obj.transform3D_2_PlaneStressInVoigt(obj.tensorVoigt);
        end
        
        
        function createRandomTensor(obj)
            obj.tensor = rand(3,3,3,3);  
            obj.MakeMajorAndMinorSymmetrization();
        end
        
        function MakeMajorAndMinorSymmetrization(obj)
            A = zeros(3,3,3,3);
            for i = 1:size(A,1)
                for j = 1:size(A,2)
                    for k = 1:size(A,3)
                        for l = 1:size(A,4)
                            Aijkl = obj.tensor(i,j,k,l);
                            Ajikl = obj.tensor(j,i,k,l);
                            Aijlk = obj.tensor(i,j,l,k);
                            Ajilk = obj.tensor(j,i,l,k);
                            Aklij = obj.tensor(k,l,i,j);
                            Alkij = obj.tensor(l,k,i,j);
                            Aklji = obj.tensor(k,l,j,i);
                            Alkji = obj.tensor(l,k,j,i);
                            A(i,j,k,l) = 1/8*( Aijkl + Ajikl + Aijlk + Ajilk + ...
                                               Aklij + Alkij + Aklji + Alkji );   
                        end
                    end
                end
            end
            obj.tensor = A;
        end
        
        function D = transform3D_2_PlaneStressInVoigt(obj,C)        
            PlaneStressVoigtTransformer = PlaneStressVoigtTensorTransformer(C);
            D = PlaneStressVoigtTransformer.tensorVoigtInPlaneStress;
        end
        
        
    end
    
    
    methods (Access = private)

        function  Cv = RepresentTensorInVoigt(obj,A)
            Cv = sym(zeros(6,6));
            for i = 1:size(A,1)
                for j = 1:size(A,2)
                    for k = 1:size(A,3)
                        for l = 1:size(A,4)
                            iv = obj.IndexTransformer.transformTensor2Voigt(i,j);
                            jv = obj.IndexTransformer.transformTensor2Voigt(k,l);
                            Cv(iv,jv) = A(i,j,k,l);
                        end
                    end
                end
            end
        end
        
    end
    
    methods (Access = protected)
 
        function vect = RepresentSymm2ndOrderTensorInVoigt(obj,A)            
            vect = sym(zeros(6,1));
            for i = 1:size(A,1)
                for j = 1:size(A,2)
                    iv = obj.transformTensorIndex2VoigtIndex(i,j);
                    vect(iv) = A(i,j);
                end
            end
        end
        
    end
    
end

