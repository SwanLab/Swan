classdef ShFunc_Chomog < Shape_Functional
    properties (Access = public)
        h_C_0
        Chomog
        tstress
        tstrain
        Chomog_Derivatives
        physicalProblem
        interpolation
        rho
        matProps
    end
    methods
        function obj=ShFunc_Chomog(settings)
            obj@Shape_Functional(settings);
            obj.physicalProblem = Physical_Problem_Micro(settings.filename);
            obj.physicalProblem.preProcess;
            obj.interpolation = Interpolation.create(settings.TOL,settings.material,settings.method);
        end
    end
    methods (Access = protected)
        function compute_Chomog_Derivatives(obj,x)
            obj.rho=obj.filter.getP0fromP1(x);
            obj.matProps=obj.interpolation.computeMatProp(obj.rho);
            %           mass=filter.Msmooth;
            %             obj.tstrain = permute(obj.tstrain,[2 3 4]);
            %             obj.tstress = permute(obj.tstrain,);
            obj.Chomog_Derivatives = zeros(obj.physicalProblem.dim.nstre,obj.physicalProblem.dim.nstre,obj.physicalProblem.geometry.ngaus,obj.physicalProblem.mesh.nelem);
            for istreChomog = 1:obj.physicalProblem.dim.nstre
                for jstreChomog = 1:obj.physicalProblem.dim.nstre
                    for igaus=1:obj.physicalProblem.geometry.ngaus
                        for istre=1:obj.physicalProblem.dim.nstre
                            for jstre = 1:obj.physicalProblem.dim.nstre
                                obj.Chomog_Derivatives(istreChomog,jstreChomog,igaus,:) = ...
                                    squeeze(obj.Chomog_Derivatives(istreChomog,jstreChomog,igaus,:)) + ...
                                    (squeeze(obj.tstrain(istreChomog,igaus,istre,:))...
                                    .*squeeze(obj.matProps.dC(istre,jstre,:,igaus))...
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
        
        function r = derivative_projection_Chomog(obj,inv_matCh,alpha,beta)
            weights = alpha*beta';
            weights_inv = inv_matCh*weights*inv_matCh;
            DtC1 = zeros(obj.physicalProblem.geometry.ngaus,obj.physicalProblem.mesh.nelem);
            DtC = zeros(obj.physicalProblem.geometry.ngaus,obj.physicalProblem.mesh.nelem);
            for igaus=1:obj.physicalProblem.geometry.ngaus
                for i=1:obj.physicalProblem.dim.nstre
                    for j=1:obj.physicalProblem.dim.nstre
                        DtC1(igaus,:) = squeeze(obj.Chomog_Derivatives(i,j,igaus,:));
                        DtC(igaus,:) = DtC(igaus,:)- weights_inv(i,j)*DtC1(igaus,:);
                    end
                end
            end
            r = DtC;
        end
        
        function computePhysicalData(obj,x)
            obj.rho=obj.filter.getP0fromP1(x);
            obj.matProps=obj.interpolation.computeMatProp(obj.rho);
            obj.physicalProblem.setMatProps(obj.matProps);
            obj.physicalProblem.computeChomog;
            obj.Chomog = obj.physicalProblem.variables.Chomog;
            obj.tstrain = obj.physicalProblem.variables.tstrain;
            obj.tstress = obj.physicalProblem.variables.tstress;
        end
    end
end







