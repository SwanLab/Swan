classdef ShFunc_Volume< Shape_Functional
    properties 
        Vfrac
    end
    methods 
        function obj=ShFunc_Volume(volumesettings)
            obj.Vfrac=volumesettings.Vfrac;
        end
        function computef(obj, x, physicalProblem, interpolation,filter)
            mass=physicalProblem.computeMass(2);
            P=filter.computePoperator(mass,physicalProblem);
            rho=filter.getP0fromP1(x,physicalProblem.mesh.coord,physicalProblem.mesh.connec,P);
            %Update phys problem
            matProps=interpolation.computeMatProp(rho);
            physicalProblem.setMatProps(matProps);
            physicalProblem.computeVariables;
            M0 = sparse(1:physicalProblem.mesh.nelem,1:physicalProblem.mesh.nelem,physicalProblem.geometry.dvolu);
            
            %compute volume
            geometric_volume = sum(mass(:));
            volume = sum(M0*rho);
            volume = volume/(geometric_volume*obj.Vfrac) - 1;
            %compute gradient
            gradient_volume = 1/(geometric_volume*obj.Vfrac);
            gradient_volume = gradient_volume*ones(size(physicalProblem.mesh.connec,1),1);
            gradient_volume = filter.getP1fromP0(gradient_volume,M0,P);
            gradient_volume = mass*gradient_volume;
            
            obj.value=volume;
            obj.gradient=gradient_volume;
        end
    end
end
