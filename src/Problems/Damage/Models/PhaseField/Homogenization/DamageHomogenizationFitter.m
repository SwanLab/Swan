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

        function fun = computeFitting(~, degPoly, phi, C)
            phi   = reshape(phi, length(phi), []);
            nStre = size(C, 1);
            fun   = cell(2,2,2,2);

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre
                            coeffs       = polyfit(phi, squeeze(C(i,j,k,l,:)), degPoly);
                            fun{i,j,k,l} = poly2sym(coeffs);

                            if isempty(symvar(fun{i,j,k,l}))
                                syms x
                                fun{i,j,k,l} = 1e-20.*x.^9;
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

% classdef DamageHomogenizationFitter < handle
% 
%     methods (Access = public, Static)
% 
%         function [fun, dfun, ddfun] = computePolynomial(~, phi, C)
%             obj = DamageHomogenizationFitter();
%             [fun, dfun, ddfun] = obj.computeSplineFitting(phi, C);
%         end
% 
%     end
% 
%     methods (Access = private)
% 
%         function [fun, dfun, ddfun] = computeSplineFitting(~, phi, C)
%             phi   = reshape(phi, length(phi), []);
%             nStre = size(C, 1);
%             fun   = cell(2,2,2,2);
%             dfun  = cell(2,2,2,2);
%             ddfun = cell(2,2,2,2);
% 
%             for i = 1:nStre
%                 for j = 1:nStre
%                     for k = 1:nStre
%                         for l = 1:nStre
%                             y = squeeze(C(i,j,k,l,:));
% 
%                             % Suavização com csaps (cubic smoothing spline)
%                             % p = 0.9 dá mais suavidade, 0.99 dá mais ajuste
%                             p = 0.99;  % ajuste este valor conforme necessário
%                             pp = csaps(phi, y, p);
% 
%                             % Derivadas usando fnder
%                             pp1 = fnder(pp, 1);
%                             pp2 = fnder(pp, 2);
% 
%                             fun{i,j,k,l}   = @(x) ppval(pp,  x);
%                             dfun{i,j,k,l}  = @(x) ppval(pp1, x);
%                             ddfun{i,j,k,l} = @(x) ppval(pp2, x);
%                         end
%                     end
%                 end
%             end
%         end
% 
%     end
% 
% end