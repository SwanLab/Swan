classdef testPrincipalDirection < testShowingError
    
     properties (Access = protected)
       tol = 1e-12;        
     end    
     
     properties (Access = protected)
        stressDim
        pdim
        nelem
        nGaus
        tensor
     end
     
    methods (Access = protected)
   
        function createTensor(obj)
            obj.tensor = rand(obj.nGaus,obj.stressDim,obj.nelem);            
        end
        
        function computeError(obj)           
            [dirS,strS] = obj.createSymbolicDirectionAndStress();
            [dirP,strP] = obj.createPrecomputedDirectionAndStress();
            error1 = norm(dirS(:) - dirP(:))/norm(dirP(:));
            error2 = norm(strS(:) - strP(:))/norm(strP(:));
            obj.error = max([error1, error2]);
        end
        
        function [dir,str] = createPrecomputedDirectionAndStress(obj)
            type = 'PRECOMPUTED';
            [dir,str] = obj.computeDirectionAndStress(type);            
        end
        
        function [dir,str] = createSymbolicDirectionAndStress(obj)
            type = 'SYMBOLIC';
            [dir,str] = obj.computeDirectionAndStress(type);            
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