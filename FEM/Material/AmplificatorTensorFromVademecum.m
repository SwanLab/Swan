classdef AmplificatorTensorFromVademecum < VariableFromVademecum
    
    
    properties (Access = protected)
        fieldName = 'Ptensor';
    end
    
    properties (Access = private)
        monomials
    end
    
    methods (Access = public)
        
        function [P,dP] = compute(obj,x)
            obj.computeParamsInfo(x);
            obj.setValuesToInterpolator(x);
            [P,dP] = obj.computeValues();
        end
        
    end
    
    methods (Access = protected)
        
        function obtainValues(obj)
            iPnorm = 1;
            obj.monomials = obj.vadVariables.monomials{iPnorm};
            var = obj.vadVariables.variables;
            mxV = obj.vadVariables.domVariables.mxV;
            myV = obj.vadVariables.domVariables.mxV;
            for imx = 1:length(mxV)
                for imy = 1:length(myV)
                    P = var{imx,imy}.(obj.fieldName);
                    Ptensor = obj.transformVector2tensor(P{1});
                    v(:,:,imx,imy) = Ptensor;
                end
            end
            obj.values = v;
        end
        
    end
    
    methods (Access = private)
        
        function [P,dP] = computeValues(obj)
            nstre = size(obj.values,1);
            P  = zeros(nstre,nstre,obj.nPoints);
            dP = zeros(nstre,nstre,obj.nPoints,obj.nParams);
            for i = 1:nstre
                for j = 1:nstre
                    pv = squeeze(obj.values(i,j,:,:));
                    [p,dp] = obj.interpolator.interpolate(pv);
                    P(i,j,:) = p;
                    dP(i,j,:,1) = dp(:,1);
                    dP(i,j,:,2) = dp(:,2);
                end
            end
        end

        function P = transformVector2tensor(obj,Pv)
            nstre = 3;
            P = zeros(nstre,nstre);
            for i = 1:nstre
                for j = 1:nstre
                    k = obj.voigt2TensorWithMonomials(i,j);
                    P(i,j) = Pv(k);
                end
            end
        end
        
    end
    
    methods (Access = private, Static)
        
        function tij = voigt2TensorWithMonomials(i,j)
            T = [6 1 2;
                 1 5 3
                 2 3 4];
            tij = T(i,j);
        end
        
    end
    
end
