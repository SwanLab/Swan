classdef testPrincipalDirection2Dnou < handle
    
    properties (Access = protected)
       tol = 1e-12;
     end
     
     properties (Access = protected)
        stressDim
        pdim
        nelem
        nGaus
        tensor
        error
     end
    
    
    properties (Access = protected)
       FileName 
    end
    
    methods (Access = public)
        
        function checkTestPassed(obj,FileName)
            obj.FileName = FileName;
            if obj.hasPassed()                
                obj.printTestPassed()
            else
                obj.printTestNotPassed()
            end
        end
        
    end

    methods (Access = public)
       
        function obj = testPrincipalDirection2Dnou()
            obj.init();
            obj.createTensor();
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj)
            obj.stressDim = 3;
            obj.pdim = '2D';
            obj.nelem = 6400; 
            obj.nGaus = 3;            
        end
        
    end
     
    methods (Access = protected)
   
        function createTensor(obj)
            obj.tensor = rand(obj.nGaus,obj.stressDim,obj.nelem);         
        end
        
        function computeError(obj)           
            [dirS,strS] = obj.computeDirectionAndStress('SYMBOLIC');
            [dirP,strP] = obj.computeDirectionAndStress('PRECOMPUTED');
            error1 = norm(dirS(:) - dirP(:))/norm(dirP(:));
            error2 = norm(strS(:) - strP(:))/norm(strP(:));
            obj.error = max([error1, error2]);
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

    methods (Access = protected)
        function printTestPassed(obj)
           cprintf('green',obj.FileName);                                    
           cprintf('green',' PASSED.');
           cprintf('black',['Error: ',num2str(obj.error),'\n']);
        end
        
        function printTestNotPassed(obj)
            cprintf('red',obj.FileName);                        
            cprintf('red',' FAILED.');
            cprintf('red',['Error: ',num2str(obj.error),'\n']);
        end
        
        function hasPassed = hasPassed(obj)
            obj.computeError()
            hasPassed = obj.error < obj.tol();
        end
    end

end