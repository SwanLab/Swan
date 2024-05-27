classdef VademecumHomogenizedVariablesLoader < handle
    
    properties (Access = private)
        variables
        structuredMesh
    end
    
    properties (Access = private)
        fileName
    end
    
    methods (Access = public)
        
        function obj = VademecumHomogenizedVariablesLoader(cParams)
            obj.init(cParams)
            obj.loadVademecum();
            obj.createStructuredMesh();            
        end

        function sM = getStructuredMesh(obj)
            sM = obj.structuredMesh;
        end

        function rho = createDensity(obj)
            rhoV = obj.loadField('volume');
            m = obj.structuredMesh.mesh;
            rho = LagrangianFunction.create(m, 1, 'P1');
            rho.fValues  = rhoV(:);
        end

        function C = createCtensor(obj)
           Cv = obj.loadField('Ctensor');
            m = obj.structuredMesh.mesh;
             for i = 1:size(Cv,3)
                 for j = 1:size(Cv,4)
                     Cij = squeeze(Cv(:,:,i,j));
                     CijF = LagrangianFunction.create(m, 1, 'P1');
                     CijF.fValues  = Cij(:);
                     C{i,j} = CijF;
                 end
             end    
        end        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.fileName    = cParams.fileName;            
        end
        
        function loadVademecum(obj)
            fName = [obj.fileName,'WithAmplificators'];
            matFile   = [fName,'.mat'];
            file2load = fullfile('Vademecums',matFile);
            v = load(file2load);
            var = v.d;            
            obj.variables = var;
        end
        
        function createStructuredMesh(obj)
            s.x = obj.variables.domVariables.mxV;
            s.y = obj.variables.domVariables.myV;
            m = StructuredMesh(s);
            obj.structuredMesh = m;
        end        
                
        function f = loadField(obj,field)
            mxV = obj.variables.domVariables.mxV;
            myV = obj.variables.domVariables.myV;
            %[] = size(obj.variables.variables{1,1}.(field))
             for imx = 1:length(mxV)
                 for imy = 1:length(myV)
                     f(imx,imy,:,:) = obj.variables.variables{imx,imy}.(field);
                 end
             end                     
        end
        
    end
    
end