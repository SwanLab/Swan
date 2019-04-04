classdef VademecumComputerForGivenVolume < handle
    
    properties (Access = protected)
      volume  
      print     
      my 
      fileName   
      prefixName            
    end
    
    
    properties (Access = private)      
      inclutionRatio
      qNorm            
    end
    
    methods (Access = public, Static)

        function obj = create(d)
            f = VademecumComputerForGivenVolumeFactory();
            obj = f.create(d);            
        end        
        
    end
                
    methods (Access = public)
        
        function compute(obj)
            obj.findInclusionLengthForCertainVolume();
            obj.obtainPrefixName();
            obj.print = true;
            obj.computeCellVariables(obj.my);            
        end
    end    
    
    methods (Access = protected)
               
        function init(obj,d)
            obj.volume          = d.volume;
            obj.inclutionRatio = d.inclutionRatio;
            obj.qNorm           = d.qNorm;
        end          
        
        function m  = computeInclusionLengthForRectangle(obj)
            m = sqrt((1-obj.volume)/obj.inclutionRatio);
        end     
        
        function computeCellVariables(obj,my)
            d = obj.computeInputForVademecumCalculator(my);
            vc = VademecumCellVariablesCalculator(d);
            vc.computeVademecumData()
            vc.saveVademecumData();
        end        
        
    end
    
    methods (Access = private)
                
        function d = computeInputForVademecumCalculator(obj,my)
            d = SettingsVademecumCellVariablesCalculator();
            d.fileName   = [obj.prefixName,obj.fileName];
            d.freeFemFileName = obj.fileName;
            d.mxMin = my*obj.inclutionRatio;
            d.mxMax = my*obj.inclutionRatio;
            d.myMin = my;
            d.myMax = my;
            d.nMx   = 1;
            d.nMy   = 1;
            d.outPutPath = [];
            d.print = obj.print;
            d.freeFemSettings.hMax = 0.02;%0.0025;
            d.freeFemSettings.qNorm = obj.qNorm;
        end
        

        
        function obtainPrefixName(obj)
            volumeStr = strrep(num2str(obj.volume),'.','');
            txiStr = strrep(num2str(obj.inclutionRatio),'.','');                                    
            obj.prefixName = ['CaseOfStudy','Rho',volumeStr,'Txi',txiStr,'QInf'];            
        end                        
        
    end
    
    methods (Access = protected, Abstract)
        findInclusionLengthForCertainVolume(obj)
    end

end