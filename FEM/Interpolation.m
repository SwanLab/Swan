classdef Interpolation < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        T
        xpoints
        isoparametric
        order
        npnod
    end
    
    methods (Static)
        function interp = create(obj)
%             class_type = class(obj);
%             if isa(obj,'Mesh')==1
%                 interp = Interpolation_geometry;
%             else
%                 
%             end
            
        switch obj
            case 'mesh'
                interp = Interpolation_geometry;
            case 'variable'
                interp = Interpolation_variable;
        end

        end
        
    end
end