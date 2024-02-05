classdef AnisotropicFromHomogenization < Material
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        vadVariables
        interpolator        
        Ctensor 
    end
    
    properties (Access = private)
        microParams
        fileName
    end
    
    methods (Access = public)
        
        function obj = AnisotropicFromHomogenization(cParams)
            obj.init(cParams)
            obj.loadVademecum();           
            obj.obtainValues();                                              
            obj.createInterpolator();
        end
        
        function C = evaluate(obj,xV)
            [C,dC] = obj.computeValues(xV);
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.microParams = cParams.microParams;
           obj.fileName = cParams.fileName;           
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
                    v(:,:,imx,imy) = var{imx,imy}.('Ctensor');
                end
            end
            obj.Ctensor = v;
        end            
        
        function createInterpolator(obj)
            sM.x = obj.vadVariables.domVariables.mxV;
            sM.y = obj.vadVariables.domVariables.myV;
            s.mesh = StructuredMesh(sM);  
            obj.interpolator = Interpolator(s);            
        end
        
        
     function [C,dCt] = computeValues(obj,xV)
            nVar = length(obj.microParams);

            mx = obj.microParams{1};
            my = obj.microParams{2} ;

            mxV = mx.evaluate(xV);
            myV = my.evaluate(xV);
            obj.interpolator.setValues(mxV(:),myV(:));            
            
            Ctens = obj.Ctensor;
            nStre = size(Ctens,1);            
            C  = zeros(nStre,nStre,mx.nDofs);
            dC = zeros(nStre,nStre,mx.nDofs,nVar);
            for i = 1:nStre
                for j = 1:nStre
                    cv = squeeze(Ctens(i,j,:,:));                    
                    [c,dc] = obj.interpolator.evaluate(cv);
                    C(i,j,:) = c;
                    for iVar = 1:nVar
                        dC(i,j,:,iVar) = dc(:,iVar);
                    end
                end
            end  
            dCt{1} = dC(:,:,:,1);
            dCt{2} = dC(:,:,:,2);
        end                  
        
    end
    
end