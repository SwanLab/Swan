classdef DensityHomogenizationFitter < handle

    methods (Access = public, Static)

        function [fun,dfun,ddfun] = computePolynomial(degPoly,rho,C)
            obj = DensityHomogenizationFitter();
            fun = obj.computeFitting(degPoly,rho,C);
            [dfun,ddfun] = obj.computeDerivative(fun);
            [fun,dfun,ddfun] = obj.convertToHandle(fun,dfun,ddfun);
        end

    end

    methods (Access = private)

        function fun = computeFitting(~,degPoly,rho,C)
            syms x

            rho = reshape(rho,length(rho),[]);

            nStre = size(C,1);

            fun = cell(2,2,2,2);

            for i=1:nStre
                for j=1:nStre
                    for k=1:nStre
                        for l=1:nStre

                            fixedPointX = [0,1];

                            Csolid = squeeze(C(i,j,k,l,end));

                            fixedPointY = [0,Csolid];

                            coeffs = polyfix(rho,squeeze(C(i,j,k,l,:)),...
                                degPoly,fixedPointX,fixedPointY);

                            fun{i,j,k,l} = poly2sym(coeffs);

                            if isempty(symvar(fun{i,j,k,l}))
                                fun{i,j,k,l} = 1e-20.*x.^degPoly;
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

                            fun{i,j,k,l} = matlabFunction(fun{i,j,k,l});

                            dfun{i,j,k,l} = matlabFunction(dfun{i,j,k,l});

                            ddfun{i,j,k,l} = matlabFunction(ddfun{i,j,k,l});

                        end
                    end
                end
            end

        end

    end

end