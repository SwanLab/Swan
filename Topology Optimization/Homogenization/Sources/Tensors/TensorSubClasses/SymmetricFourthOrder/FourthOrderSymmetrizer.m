classdef FourthOrderSymmetrizer < handle

    properties
    end

    methods

        function obj = FourthOrderSymmetrizer()
        end
    end

    methods (Access = public)

        function isSymmetric = isSymmetric(obj,Tensor)
            Tsym = obj.symmetrize(Tensor);
            isSymmetric = norm(Tsym(:) - Tensor(:)) < 1e-12;
        end

    end

    methods (Access = public, Static)

        function A = symmetrize(Tensor)
            T = Tensor;
            A = zeros(3,3,3,3);
            for i = 1:size(A,1)
                for j = 1:size(A,2)
                    for k = 1:size(A,3)
                        for l = 1:size(A,4)
                            Aijkl = T(i,j,k,l);
                            Ajikl = T(j,i,k,l);
                            Aijlk = T(i,j,l,k);
                            Ajilk = T(j,i,l,k);
                            Aklij = T(k,l,i,j);
                            Alkij = T(l,k,i,j);
                            Aklji = T(k,l,j,i);
                            Alkji = T(l,k,j,i);
                            A(i,j,k,l) = 1/8*( Aijkl + Ajikl + Aijlk + Ajilk + ...
                                Aklij + Alkij + Aklji + Alkji );
                        end
                    end
                end
            end

        end

    end



end

