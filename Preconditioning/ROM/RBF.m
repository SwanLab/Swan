classdef RBF < handle
    % RADIAL BASIS FUNCTION
    properties (Access = private)
        Xdata
        Ydata
        centers
        sigma
    end

    properties (Access = private)
        f
    end

    methods (Access = public)

        function obj = RBF(Ydata,Xdata,centers)
            obj.init(Ydata,Xdata,centers)
            obj.constructFunction();
        end

        function y = evaluate(obj,r)
            y = obj.f(r);
        end

    end

    methods (Access = private)

        function init(obj,Ydata,Xdata,centers)
            obj.Ydata   = Ydata;
            obj.Xdata   = Xdata;
            obj.centers = centers;
            center_spacing = mean(diff(centers));
            nbasis = size(Ydata,2);
            obj.sigma = 0.8*center_spacing*sqrt(nbasis);
        end

        function constructFunction(obj)
            R = (obj.Xdata - obj.centers').^2;
            Phi = exp(-R / (2 * obj.sigma^2));
            w = Phi\ obj.Ydata;
            obj.f = @(x)  w'*exp(-((x - obj.centers').^2) / (2 * obj.sigma^2));
%             x_plot = linspace(min(obj.Xdata), max(obj.Xdata), 500);
%             for i = 1:length(obj.centers)
%                 Phi_plot(:, i) = exp(-((x_plot - obj.centers(i)).^2) / (2 * obj.sigma^2));
%             end
%             figure;
%             plot( Phi_plot, 'LineWidth', 1.5);
        end

    end

end