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
            obj.sigma = center_spacing / sqrt(2);
        end

        function constructFunction(obj)
            R = (obj.Xdata - obj.centers').^2;
            Phi = exp(-R / (2 * obj.sigma^2));
            w = Phi\ obj.Ydata;
            obj.f = @(x)  w'*exp(-((x - obj.centers').^2) / (2 * obj.sigma^2));
        end

    end

end