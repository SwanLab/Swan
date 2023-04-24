classdef DesignVariable < handle
    
    properties (Access = public)
        thetavec
    end

    properties (Access = private)
        neuronsPerLayer
        nLayers
    end
    
    methods (Access = public)

        function obj = DesignVariable(s)
            obj.init(s);
        end

    end

    methods (Access = private)

        function init(obj,s)
           obj.thetavec = s.initValue;
       end  
    end
end

