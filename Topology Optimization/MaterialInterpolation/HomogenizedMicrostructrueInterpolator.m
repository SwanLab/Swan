classdef HomogenizedMicrostructrueInterpolator < handle

    properties (Access = private)
       Ctensor 
    end
    
    properties (Access = private)
        fileName
        vadVariables
        interpolator
        
    end
    
    methods (Access = public)
        
        function obj = HomogenizedMicrostructrueInterpolator(cParams)
            obj.init(cParams)
            obj.loadVademecum();           
            obj.obtainValues();                      
            obj.createInterpolator();
        end

        function [Cm,dCm] = computeConsitutiveTensor(obj,x)
 %           obj.nVar = length(x);

%            Cref  = c;
%            dCref = permute(dc,[1 2 4 3 5]);

            [C,dC] = obj.computeValues(x);   

         %   xR = obj.reshapeDesignVariable(x);
            Cm = obj.createMaterial(C);
            for iVar = 1:length(x)
               dCi = dC(:,:,:,iVar);
               dCm{iVar} = obj.createMaterial(dCi);
            end
        end
        
    end
 
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.fileName = cParams.fileName;
           obj.microParams = cParams.microParams;
        end

        function createInterpolator(obj)
            sM.x = obj.vadVariables.domVariables.mxV;
            sM.y = obj.vadVariables.domVariables.myV;
            s.mesh = StructuredMesh(sM);  
            obj.interpolator = Interpolator(s);            
        end

        function [C,dCt] = computeValues(obj,x)
            nVar = length(x);
            mx = x{1};
            my = x{2}; 
         %   nValues = mx.nn
            obj.interpolator.setValues(mx.fValues,my.fValues);            
            
            Ctens = obj.Ctensor;
            nStre = size(Ctens,1);            
            C  = zeros(nStre,nStre,mx.nDofs);
            dC = zeros(nStre,nStre,mx.nDofs,nVar);
            for i = 1:nStre
                for j = 1:nStre
                    cv = squeeze(Ctens(i,j,:,:));                    
                    [c,dc] = obj.interpolator.interpolate(cv);
                    C(i,j,:) = c;
                    for iVar = 1:nVar
                        dC(i,j,:,iVar) = dc(:,iVar);
                    end
                end
            end  
            dCt{1} = dC(:,:,:,1);
            dCt{2} = dC(:,:,:,2);
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
        
        function loadVademecum(obj)
            fName = [obj.fileName,'WithAmplificators'];
            matFile   = [fName,'.mat'];
            file2load = fullfile('Vademecums',matFile);
            v = load(file2load);
            obj.vadVariables = v.d;          
        end

        function m = createMaterial(obj,C)
            s.type    = 'ANISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = 2;
            s.constiutiveTensor = C;            
            s.bulk    = kappa;
            m = Material.create(s);   
        end
        
    end
    
    methods (Access = private, Static)
        
        function xR = reshapeDesignVariable(x)
            xV = x.value;
            nV = x.nVariables;
            nx = length(xV)/nV;
            xR = cell(nV,1);
            for ivar = 1:nV
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xR{ivar} = xV(i0:iF);
            end
        end
        
    end    
    
end