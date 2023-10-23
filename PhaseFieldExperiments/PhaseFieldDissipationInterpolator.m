classdef PhaseFieldDissipationInterpolator < handle

    properties (Access = public)

    end

    properties (Access = private)
        alpha
        dalpha
        ddalpha
    end

    properties (Access = private)
       constitutiveProperties        
    end

    methods (Access = public)

        function obj = PhaseFieldDissipationInterpolator()   
            obj.createDissipationInterpolation()
        end

        function alpha = computeAlphaProp(obj,phi)
            alpha = obj.alpha(phi);
        end

        function dalpha = computeDAlphaProp(obj,phi)
            dalpha = obj.dalpha(phi);
        end

        function ddalpha = computeDDAlphaProp(obj,phi)
            ddalpha = obj.ddalpha(phi);
        end

    end

    methods (Access = private)
        function createDissipationInterpolation(obj)
            p = 2;

            obj.alpha    = @(phi) phi.^p;
            obj.dalpha   = @(phi) p*(phi).^(p-1);
            obj.ddalpha  = @(phi) p*(p-1)*(phi).^(p-2);
        end

    end

end