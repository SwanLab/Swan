classdef PrincipalDirectionTest < handle
     
     properties (Access = private)
        stressDim
        pdim
        mesh
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
            obj.mesh     = cParams.mesh;
            obj.nGaus     = cParams.nGaus;
            obj.stressDim = cParams.stressDim;
            obj.createTensor();
        end
        
    end

    methods (Access = public)

        function error = computeError(obj)
            [dirS,strS] = obj.computeDirectionAndStress('SYMBOLIC');
            [dirP,strP] = obj.computeDirectionAndStress('PRECOMPUTED');
            for iCell = 1:numel(dirS)
                errorDir(iCell) = Norm(dirS{iCell}-dirP{iCell},'L2')/Norm(dirP{iCell},'L2');
            end
            errorStr = Norm(strS -strP,'L2')/Norm(strP,'L2');            
            error = max([norm(errorDir),errorStr]);
        end

    end

    methods (Access = protected)

        function createTensor(obj)
            s.fValues = rand(obj.mesh.nelem,obj.stressDim);
            s.mesh    = obj.mesh;
            s.order   = 'P0';
            obj.tensor = LagrangianFunction(s);
        end            
        
        
        function [dir,str] = computeDirectionAndStress(obj,type)
            s.eigenValueComputer.type = type;
            s.type = obj.pdim;
            p = PrincipalDirectionComputer.create(s);
            [dir,str] = p.compute(obj.tensor);
        end

    end

end