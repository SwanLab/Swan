classdef Triangle_Constant < Isoparametric
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = Triangle_Constant
            obj = obj@Isoparametric;
            obj.type = 'TRIANGLE';
            obj.ndime = 2;
            obj.nnode = 1;
%             obj.ngaus = 1;
%             obj.weigp = 1/2;
%             obj.posgp = [1/3;1/3];
            obj.pos_nodes = [1/3 1/3]; %[0 0; 1 0; 0 1];
            
            obj.shape = @(s,t) {1};
            
            obj.deriv = @(s,t) {0;0};

        end
        
    end
end
