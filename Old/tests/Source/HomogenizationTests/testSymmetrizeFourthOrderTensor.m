classdef testSymmetrizeFourthOrderTensor < handle
    
    properties (Access = protected)
           Cani
           Csym
           energy
           nonEnergy
    end
    
    methods (Access = public)
        
        function obj = testSymmetrizeFourthOrderTensor()
            obj.computeAnisotropicFourthOrderTensor();
            obj.computeAplicationWithRandomSymetricTensor()
        end

        function hasPassed = hasPassed(obj)
            meanEn = mean(obj.energy(1:4));
            desV = obj.energy(1:4) - meanEn;
            PreservedBySymmetry = norm(desV) < 1e-12;
            NotPreservedByNoSymmetry = abs(meanEn - obj.nonEnergy) > 1e-12;
            hasPassed = PreservedBySymmetry && NotPreservedByNoSymmetry ;
        end

    end

    methods (Access = private)

        function computeAnisotropicFourthOrderTensor(obj)
            obj.Cani = SymmetricFourthOrder3DTensor();
            obj.Cani.createRandomTensor();
        end
        
        function computeAplicationWithRandomSymetricTensor(obj)
            strain = Strain3DTensor;
            strain.createRandomTensor();
            txi = strain.getValue();
            A = obj.Cani.getValue();
            e = zeros(4,1);
            nE = 0;
            for i = 1:size(A,1)
                for j = 1:size(A,2)
                    for k = 1:size(A,3)
                        for l = 1:size(A,4)
                            e(1) = e(1) + txi(i,j)*A(i,j,k,l)*txi(k,l);
                            e(2) = e(2) + txi(j,i)*A(i,j,k,l)*txi(k,l);
                            e(3) = e(3) + txi(i,j)*A(i,j,k,l)*txi(l,k);
                            e(4) = e(4) + txi(k,l)*A(i,j,k,l)*txi(i,j);
                            nE = nE + txi(i,k)*A(i,j,k,l)*txi(j,l);
                        end
                    end
                end
            end
            obj.energy = e;
            obj.nonEnergy = nE;
        end

    end

end