classdef AnisotropicFromHomogenization < Material
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        vadVariables
        Ctensor 
        sMesh
    end
    
    properties (Access = private)
        microParams
        fileName
    end
    
    methods (Access = public)
        
        function obj = AnisotropicFromHomogenization(cParams)
            obj.init(cParams)
            obj.loadVademecum();    
            obj.createStructuredMesh();
            obj.obtainValues();      
        end
        
        function C = evaluate(obj,xV)
            C = obj.computeValues(xV);            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.microParams = cParams.microParams;
           obj.fileName    = cParams.fileName;           
        end
        
        function loadVademecum(obj)
            fName = [obj.fileName,'WithAmplificators'];
            matFile   = [fName,'.mat'];
            file2load = fullfile('Vademecums',matFile);
            v = load(file2load);
            obj.vadVariables = v.d;          
        end
        
        function obtainValues(obj)
            var = obj.vadVariables.variables;
            mxV = obj.vadVariables.domVariables.mxV;
            myV = obj.vadVariables.domVariables.mxV;
            for imx = 1:length(mxV)
                for imy = 1:length(myV)
                    C(:,:,imx,imy) = var{imx,imy}.('Ctensor');
                end
            end   
            
            for i = 1:size(C,1)
                for j = 1:size(C,2)
                    Cij = squeeze(C(i,j,:,:));
                    s.fValues = Cij(:);
                    s.mesh    = obj.sMesh.mesh;
                    s.ndim    = 1;
                    CijF = LagrangianFunction.create(obj.sMesh.mesh, 1, 'P1');
                    obj.Ctensor{i,j} = CijF;
                end
            end
            
        end         

        function createStructuredMesh(obj)
            s.x = obj.vadVariables.domVariables.mxV;
            s.y = obj.vadVariables.domVariables.myV;
            m = StructuredMesh(s); 
            obj.sMesh = m;
        end     

        function [mL,cells] = obtainLocalCoord(obj,xV)
            mx = obj.microParams{1};
            my = obj.microParams{2};
            mxG = mx.evaluate(xV);
            myG = my.evaluate(xV);
            mG(:,1) = mxG(:);
            mG(:,2) = myG(:);
            [mL,cells] = obj.sMesh.obtainLocalFromGlobalCoord(mG);
        end
        
        function C = computeValues(obj,xV)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nStre = size(obj.Ctensor,1); 
            nDofs = size(mL,2);
            C  = zeros(nStre,nStre,nDofs);
            for i = 1:nStre
                for j = 1:nStre
                    Cij(1,1,:) = obj.Ctensor{i,j}.sample(mL,cells);  
                    C(i,j,:)   = Cij;
                end
            end
        end                  
        
    end
    
end