classdef DamageHomogenizationFitter < handle

    methods (Access = public, Static)

        function [fun,dfun,ddfun] = computePolynomial(degPoly,phi,C)
            obj = DamageHomogenizationFitter();
            fun = obj.computeFitting(degPoly,phi,C);
            [dfun,ddfun] = obj.computeDerivative(fun);
            [fun,dfun,ddfun] = obj.convertToHandle(fun,dfun,ddfun);
        end

    end

    methods (Access = private)

        function fun = computeFitting(~,degPoly,phi,C)
            syms x
            phi = reshape(phi,length(phi),[]);

            nStre = size(C,1);
            fun   = cell(3,3);
            for i=1:nStre
                for j=1:nStre
                    fixedPointX = [0,1];
                    fixedPointY = [squeeze(C(i,j,1)),0];
                    coeffs = polyfix(phi,squeeze(C(i,j,:)),degPoly,fixedPointX,fixedPointY);
                    fun{i,j} = poly2sym(coeffs);
                    if isempty(symvar(fun{i,j}))
                        fun{i,j} = 1e-20.*x.^9;
                    end
                end
            end
        end

        function [dfun,ddfun] = computeDerivative(~,fun)
            nStre = size(fun,1);
            dfun  = cell(3,3);
            ddfun = cell(3,3);
            for i=1:nStre
                for j=1:nStre
                    dfun{i,j} = diff(fun{i,j});
                    ddfun{i,j} = diff(dfun{i,j});
                end
            end
        end

        function [fun,dfun,ddfun] = convertToHandle(~,fun,dfun,ddfun)
            nStre = size(fun,1);
            for i=1:nStre
                for j=1:nStre
                        fun{i,j}   = matlabFunction(fun{i,j});
                        dfun{i,j}  = matlabFunction(dfun{i,j});
                        ddfun{i,j} = matlabFunction(ddfun{i,j});
                end
            end
        end
    end

end