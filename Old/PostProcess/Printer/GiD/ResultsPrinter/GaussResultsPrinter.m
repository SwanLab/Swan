classdef GaussResultsPrinter < handle
    
    properties (Abstract, Access = protected)
        dataBase
    end
    
    methods (Access = public)
        
        function obj = GaussResultsPrinter(d,dh)

        end
        
        function dB = getDataBase(obj)
            dB = obj.dataBase;
        end
         
    end
    
    methods (Access = protected)

        function createHeadPrinterDataBase(obj,d,dh)
            dI.fileID = dh.fileID;
            dI.etype  = dh.etype;
            dI.ndim   = dh.ndim;
            dI.gaussDescriptor = 'Guass up?';
            dI.posgp = d.quad.posgp';
            dI.ngaus = d.quad.ngaus;
            obj.headPrinterDataBase = dI;
        end
        
%         function storeQuadInfo(obj,d)
%             obj.ngaus = d.quad.ngaus;
%             obj.posgp = d.quad.posgp';
%         end
        
    end
    
%     methods (Access = protected, Abstract)
%         createScalarDataBase(obj)
%         createVectorDataBase(obj)
%     end
    
end