classdef ShFunc_Chomog < Shape_Functional
    properties (Access = protected)
        h_C_0;
        Chomog;
        tstress;
        tstrain;
        Chomog_Derivatives;
    end
    methods
        function obj=ShFunc_Chomog(settings)
            obj@Shape_Functional(settings);
        end
    end
    methods (Access = protected)
        function compute_Chomog_Derivatives(obj,nstre,nelem,ngaus,x,interpolation,filter)
            rho=filter.getP0fromP1(x);
            matProps=interpolation.computeMatProp(rho);
            mass=filter.Msmooth;
            
            obj.Chomog_Derivatives = zeros(nstre,nstre,ngaus,nelem);
            for istreChomog = 1:nstre
                for jstreChomog = 1:nstre
                    for igaus=1:ngaus
                        for istre=1:nstre
                            for jstre = 1:nstre
                                obj.Chomog_Derivatives(istreChomog,jstreChomog,igaus,:) = ...
                                    squeeze(obj.Chomog_Derivatives(istreChomog,jstreChomog,igaus,:)) + ...
                                    (squeeze(obj.tstrain(istreChomog,igaus,istre,:))...%(squeeze(-obj.tstrain(istreChomog,igaus,istre,:))...
                                    .*squeeze(matProps.dC(istre,jstre,:,igaus))...
                                    .*squeeze(obj.tstrain(jstreChomog,igaus,jstre,:)));
                            end
                        end
                    end
%                     C_D = filter.getP1fromP0(squeeze(obj.Chomog_Derivatives(istreChomog,jstreChomog,:,:)));
%                     obj.Chomog_Derivatives(istreChomog,jstreChomog,:,:) = mass*C_D;
                end
            end
            
        end
        
        function r = projection_Chomog(obj,inv_matCh,alpha,beta)
            weights = alpha*beta';
            r = sum(sum(weights.*inv_matCh));
        end
        
        function r = derivative_projection_Chomog(obj,inv_matCh,alpha,beta,Chomog_Derivatives,nelem,ngaus,nstre)
            weights = alpha*beta';
            weights_inv = inv_matCh*weights*inv_matCh;
            DtC1 = zeros(ngaus,nelem);
            DtC = zeros(ngaus,nelem);
            for igaus=1:ngaus
                for i=1:nstre
                    for j=1:nstre
                        DtC1(igaus,:) = squeeze(Chomog_Derivatives(i,j,igaus,:));
                        DtC(igaus,:) = DtC(igaus,:)- weights_inv(i,j)*DtC1(igaus,:);
                    end
                end
            end
            r = DtC;
        end
    end
end







