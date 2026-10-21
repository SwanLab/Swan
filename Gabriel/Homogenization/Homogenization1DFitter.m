classdef Homogenization1DFitter < handle

    methods (Access = public, Static)

        function [fun,dfun,ddfun] = computePolynomial(degPoly,xValues,C)

            obj = Homogenization1DFitter();

            obj.validateInput(degPoly,xValues,C);

            fun = obj.computeFitting( ...
                degPoly,xValues,C);

            [dfun,ddfun] = ...
                obj.computeDerivative(fun);

            [fun,dfun,ddfun] = ...
                obj.convertToHandle( ...
                fun,dfun,ddfun);

        end

    end

    methods (Access = private)

        function validateInput(~,degPoly,xValues,C)

            xValues = xValues(:);

            nPoints = numel(xValues);

            if degPoly >= nPoints
                error( ...
                    'Homogenization1DFitter:InvalidPolynomialDegree', ...
                    'Polynomial degree must be smaller than the number of data points.');
            end

            if size(C,5) ~= nPoints
                error( ...
                    'Homogenization1DFitter:InvalidNumberOfSamples', ...
                    'The fifth dimension of C must match the number of parameter values.');
            end

        end

        function fun = computeFitting(~,degPoly,xValues,C)

            syms x

            xValues = xValues(:);

            nStre = size(C,1);

            fun = cell(2,2,2,2);

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre

                            values = ...
                                squeeze(C(i,j,k,l,:));

                            values = ...
                                values(:);

                            coeffs = ...
                                polyfit( ...
                                xValues, ...
                                values, ...
                                degPoly);

                            fun{i,j,k,l} = ...
                                poly2sym(coeffs,x);

                        end
                    end
                end
            end

        end

        function [dfun,ddfun] = ...
                computeDerivative(~,fun)

            syms x

            nStre = size(fun,1);

            dfun  = cell(2,2,2,2);
            ddfun = cell(2,2,2,2);

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre

                            dfun{i,j,k,l} = ...
                                diff( ...
                                fun{i,j,k,l},x);

                            ddfun{i,j,k,l} = ...
                                diff( ...
                                dfun{i,j,k,l},x);

                        end
                    end
                end
            end

        end

        function [fun,dfun,ddfun] = ...
                convertToHandle(~,fun,dfun,ddfun)

            syms x

            nStre = size(fun,1);

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre

                            fun{i,j,k,l} = ...
                                matlabFunction( ...
                                fun{i,j,k,l}, ...
                                'Vars',x);

                            dfun{i,j,k,l} = ...
                                matlabFunction( ...
                                dfun{i,j,k,l}, ...
                                'Vars',x);

                            ddfun{i,j,k,l} = ...
                                matlabFunction( ...
                                ddfun{i,j,k,l}, ...
                                'Vars',x);

                        end
                    end
                end
            end

        end

    end

end