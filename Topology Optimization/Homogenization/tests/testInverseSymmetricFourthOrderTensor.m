classdef testInverseSymmetricFourthOrderTensor < test
    
    properties (Access = private)
        Dim
        Tensor
        InvTensor
        Identity
        ExpectedIdentity
    end
       
    
    methods (Access = public)
        
        function obj = testInverseSymmetricFourthOrderTensor() 
            obj.init()
            obj.createRandomFourthOrderTensor()
            obj.invertFourthOrderTensor()
            obj.createIdentityFourthOrderTensor()
            obj.obtainExpectedIdentityFourthOrderTensor()
        end
      
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.Dim = 3;
        end
        
        function createRandomFourthOrderTensor(obj)
            d = obj.Dim;
            obj.Tensor = rand(d,d,d,d);
            obj.Tensor = FourthOrderSymmetrizer.symmetrize(obj.Tensor);
        end
        
        function invertFourthOrderTensor(obj)
            A = FourthOrderTensor();
            A.setValue(obj.Tensor);
            Ainv = Inverter.invert(A);
            obj.InvTensor = Ainv;
        end
        
        function createIdentityFourthOrderTensor(obj)
            d = obj.Dim;
            Id = zeros(d,d,d,d);
            I = eye(d);
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l = 1:d
                            Id(i,j,k,l) = 0.5*(I(i,k)*I(j,l)+I(i,l)*I(j,k));
                        end
                    end
                end
            end
            obj.Identity = Id;
        end
        
        function obtainExpectedIdentityFourthOrderTensor(obj)
            d = obj.Dim;
            Id = zeros(d,d,d,d);
            A = obj.Tensor;
            Ainv = obj.InvTensor;
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l = 1:d
                            for m = 1:d
                                for n = 1:d
                                    Id(i,j,k,l) = Id(i,j,k,l) + A(i,j,m,n)*Ainv(m,n,k,l);
                                end
                            end
                            
                        end
                    end
                end
            end
            obj.ExpectedIdentity = Id;
        end
        
    end
    
    
    methods (Access = protected)
                
        function hasPassed = hasPassed(obj)
            I = obj.Identity;
            Ie = obj.ExpectedIdentity;            
            hasPassed = norm(I(:)-Ie(:))/norm(I(:)) < 1e-13;
        end
    end
    
end
