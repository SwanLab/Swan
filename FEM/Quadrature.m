classdef Quadrature < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        posgp
        weigp
        ngaus
    end    
    methods
        function obj=Quadrature(order)
            obj.computeQuadrature(order)            
        end
        computeQuadrature(oborder)
    end
end
