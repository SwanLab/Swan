classdef VademecumDataLoader < handle
    
   properties (Access = public)
      cellDataSameSign 
      cellDataDifferentSign
   end
   
   properties (Access = private)
        vademecumPath
        fileNameSameSign
        fileNameDiffSign      
        nMx 
        nMy
        nPhi
   end    
    
   methods (Access = public)
       
       function obj = VademecumDataLoader(cParams)
            obj.init(cParams);
            obj.loadVademecumDataBase();
       end
       
   end

   methods (Access = private)
       
       function init(obj,cParams)
          obj.nMx = cParams.nMx;
          obj.nMy = cParams.nMy;
          obj.nPhi = cParams.nPhi;
          obj.vademecumPath = '/media/alex/My Passport/Vademecum/';
          obj.fileNameSameSign = 'OptimalSuperEllipseSameStressSign';
          obj.fileNameDiffSign = 'OptimalSuperEllipseDifferentStressSign';                    
       end    
       
        function loadVademecumDataBase(obj)
            cS = cell(obj.nMx,1);
            cD = cell(obj.nMx,1);
            for iMx = 1:obj.nMx
                [cS{iMx},cD{iMx}] = obj.loadBothVademecum(iMx);
            end
            obj.cellDataSameSign      = cS;
            obj.cellDataDifferentSign = cD;
        end
        
        function [cS,cD] = loadBothVademecum(obj,iMx)
            cS = obj.loadSameSignVademecum(iMx);
            cD = obj.loadDifferentSignVademecum(iMx);
        end
        
        function d = loadSameSignVademecum(obj,iMx)
            fName = obj.fileNameSameSign;
            d = obj.loadVademecum(fName,iMx);
        end
        
        function d = loadDifferentSignVademecum(obj,iMx)
            fName = obj.fileNameDiffSign;
            d = obj.loadVademecum(fName,iMx);
        end       
        
     function cS = loadVademecum(obj,fileName,iMx)
            fName = [obj.vademecumPath,fileName,num2str(iMx),'.mat'];
            d = load(fName);
            cS = cell(obj.nMx,obj.nMy,obj.nPhi);
            for iMy = 1:obj.nMy
                for iphi = 1:obj.nPhi
                    cS{iMx,iMy,iphi}.mesh = d.c{iMx,iMy,iphi}.mesh;
                    cS{iMx,iMy,iphi}.rho = d.c{iMx,iMy,iphi}.rho;
                    cS{iMx,iMy,iphi}.q   = d.c{iMx,iMy,iphi}.q;
                    cS{iMx,iMy,iphi}.xi  = d.c{iMx,iMy,iphi}.xi;
                    cS{iMx,iMy,iphi}.phi = d.c{iMx,iMy,iphi}.phi;
                end
            end
        end                
       
   end
    
    
end