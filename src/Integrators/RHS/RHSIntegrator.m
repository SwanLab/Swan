classdef RHSIntegrator < handle

    properties (Access = protected)
        mesh
        quadrature
        quadratureOrder
        test
    end
   

    methods (Access = public, Static)
        
        function obj = create(s)
            f = RHSIntegratorFactory();
            obj = f.create(s);
        end

    end
        
end