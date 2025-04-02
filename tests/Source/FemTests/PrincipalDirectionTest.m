classdef PrincipalDirectionTest < handle
     
     properties (Access = protected)
        stressDim
        pdim
        nelem
        nGaus
        tensor
     end
    
    
    properties (Access = protected)
       FileName
    end
    
    methods (Access = public)
       
        function obj = PrincipalDirectionTest(cParams)
            obj.init(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj, cParams)
            obj.pdim      = cParams.pdim;
            obj.nelem     = cParams.nelem;
            obj.nGaus     = cParams.nGaus;
            obj.stressDim = cParams.stressDim;
            obj.createTensor();
        end
        
    end
    
    methods (Access = public)

        function error = computeError(obj)
            [dirS,strS] = obj.computeDirectionAndStress('SYMBOLIC');
            [dirP,strP] = obj.computeDirectionAndStress('PRECOMPUTED');
            errorDir = norm(dirS(:) - dirP(:))/norm(dirP(:));
            errorStr = norm(strS(:) - strP(:))/norm(strP(:));
            error = max([errorDir, errorStr]);
        end

    end

    methods (Access = protected)
   
        function createTensor(obj)
            obj.tensor = rand(obj.nGaus,obj.stressDim,obj.nelem);
        end
        
        function [dir,str] = computeDirectionAndStress(obj,type)
            s.eigenValueComputer.type = type;
            s.type = obj.pdim;
            p = PrincipalDirectionComputer.create(s);
            p.compute(obj.tensor);
            dir = p.direction;
            str = p.principalStress;
        end

    end

end