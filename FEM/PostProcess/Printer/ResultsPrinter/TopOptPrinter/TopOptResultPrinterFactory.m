classdef TopOptResultPrinterFactory < handle
        
   methods (Access = public, Static)
       
       function p = create(d,dT,hasGaussData)
           
           if hasGaussData
               p = TopOptResultPrinterWithGauss(d,dT);
           else
               p = TopOptResultPrinterWithNoGauss(d,dT);
           end                      
       end       
       
   end

    
end