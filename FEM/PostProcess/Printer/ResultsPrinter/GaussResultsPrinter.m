classdef GaussResultsPrinter < handle
    
    properties (Access = protected)
        ngaus
        posgp
        gaussDescriptor = 'Guass up?';
        hasGaussData = true
    end
    
    
    methods (Access = public)
        
        function obj = GaussResultsPrinter()
            obj.hasGaussData = true;
        end
        
    end
    
    methods (Access = protected)
        function d = createScalarGaussDataBase(obj,varargin)
            d = obj.createScalarDataBase(varargin{:});
            d = obj.addGaussDescriptor(d);
        end
        
        function d = createVectorGaussDataBase(obj,varargin)
            d = obj.createVectorDataBase(varargin{:});
            d = obj.addGaussDescriptor(d);
        end
        
        function storeQuadInfo(obj,d)
            obj.ngaus = d.quad.ngaus;
            obj.posgp = d.quad.posgp';
        end
        
    end
    
    methods (Access = protected, Abstract)
        createScalarDataBase(obj)
        createVectorDataBase(obj)
    end
    
    methods (Access = private)
        
        function d = addGaussDescriptor(obj,d)
            d.gaussDescriptor = obj.gaussDescriptor;
        end
        
    end
    
    
end