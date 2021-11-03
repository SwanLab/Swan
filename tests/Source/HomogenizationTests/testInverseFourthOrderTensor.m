classdef testInverseFourthOrderTensor < handle
    
    properties (Access = public)
        tol = 1e-12;
    end

    properties (Access = protected)
        tensor
    end
    
    properties (Access = private)
        invTensor
        identity
        expectedIdentity
    end

    methods (Access = public)
        
        function obj = testInverseFourthOrderTensor() 
            obj.createRandomFourthOrderTensor()
            obj.invertFourthOrderTensor()
            obj.createIdentityFourthOrderTensor()
            obj.obtainExpectedIdentityFourthOrderTensor()
        end

        function error = computeError(obj)
            I = obj.identity;
            Ie = obj.expectedIdentity;            
            error = norm(I(:)-Ie(:))/norm(I(:));
        end

    end
    
    methods (Access = private)

        function invertFourthOrderTensor(obj)
            A = obj.tensor;
            Ainv = Inverter.invert(A);
            obj.invTensor = Ainv;
        end
        
        function createIdentityFourthOrderTensor(obj)
            d = obj.tensor.getDimension();
            Id = zeros(d,d,d,d);
            I = eye(d);
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l = 1:d
                            Id(i,j,k,l) = obj.computeIdentityTensor(I,i,j,k,l);
                        end
                    end
                end
            end
            obj.identity = Id;
        end
        
        function obtainExpectedIdentityFourthOrderTensor(obj)
            d = obj.tensor.getDimension();
            Id = zeros(d,d,d,d);
            A = obj.tensor.getValue();
            Ainv = obj.invTensor.getValue();
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
            obj.expectedIdentity = Id;
        end
        
    end

    methods (Abstract, Access = protected, Static)
        computeIdentityTensor(obj)
        createRandomFourthOrderTensor(obj)
    end
    
end