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

        function initDeriv = computeInitialDerivative(~)
            Gc=5e-3; l0=0.1; E=210; cw = 3/8; sigCrit =1;
            initDeriv = -2*cw*(Gc/l0)*E*(1/sigCrit)^2;
        end

        function fun = computeFitting(~,degPoly,phi,C)
            syms x
            phi = reshape(phi,length(phi),[]);

            nStre = size(C,1);
            fun   = cell(2,2,2,2);
            for i=1:nStre
                for j=1:nStre
                    for k=1:nStre
                        for l=1:nStre
                            if  (i==1 && j==1 && k~=l) || (i==2 && j==2 && k~=l) || (i~=j && k==l)
                                fun{i,j,k,l} = 1e-20.*x.^9;
                            else
                                fixedPointX = [0,1];
                                if i==1 && j==1 && k==1 && l==1
                                    fixedPointY = [squeeze(C(i,j,k,l,1)),squeeze(C(i,j,k,l,end))];
                                else
                                    fixedPointY = [squeeze(C(i,j,k,l,1)),0];
                                end
                                coeffs = polyfix(phi,squeeze(C(i,j,k,l,:)),degPoly,fixedPointX,fixedPointY);
                                fun{i,j,k,l} = poly2sym(coeffs);
                            end
                        end
                    end
                end
            end
        end

        function [dfun,ddfun] = computeDerivative(~,fun)
            nStre = size(fun,1);
            dfun  = cell(2,2,2,2);
            ddfun = cell(2,2,2,2);
            for i=1:nStre
                for j=1:nStre
                    for k=1:nStre
                        for l=1:nStre
                            dfun{i,j,k,l} = diff(fun{i,j,k,l});
                            ddfun{i,j,k,l} = diff(dfun{i,j,k,l});
                        end
                    end
                end
            end
        end

        function [fun,dfun,ddfun] = convertToHandle(~,fun,dfun,ddfun)
            nStre = size(fun,1);
            for i=1:nStre
                for j=1:nStre
                    for k=1:nStre
                        for l=1:nStre
                            fun{i,j,k,l}   = matlabFunction(fun{i,j,k,l});
                            dfun{i,j,k,l}  = matlabFunction(dfun{i,j,k,l});
                            ddfun{i,j,k,l} = matlabFunction(ddfun{i,j,k,l});
                        end
                    end
                end
            end
        end
    end

end