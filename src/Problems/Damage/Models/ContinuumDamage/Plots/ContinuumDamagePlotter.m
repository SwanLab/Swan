classdef ContinuumDamagePlotter < handle
    
    properties (Access = private)
        damage
        damageField
        displacement
        reaction
        qMax
        rMax
        energy
        iter
    end
    
    methods (Access = public)
        
        function obj = ContinuumDamagePlotter(cParams)
            obj.init(cParams)
            obj.plotForceDisplacement()
            obj.plotDamage()
            obj.meshDamage()
            obj.plotQR()
            obj.plotEnergy()
            obj.plotIterations()
        end
        
    end
        
    methods (Access = private)
        
        function init(obj,cParams)
            obj.damage = cParams.damage.maxValue;
            obj.damageField = cParams.damage.field{end};
            obj.displacement = cParams.displacement.value;
            obj.reaction = cParams.reaction;
            obj.qMax = cParams.qMaxValue;
            obj.rMax = cParams.rMaxValue;
            obj.energy = cParams.energy;
            obj.iter = cParams.iter;
        end

        function plotForceDisplacement(obj)
            figure()
            plot(obj.displacement,obj.reaction)
            title('Force-displacement diagram (Reaction force)')
            grid on
            xlabel('Displacement [mm]')
            ylabel('Force [kN]')
        end
        
        function plotDamage(obj)
            figure()
            plot(obj.displacement,obj.damage)
            title('Damage-displacement diagram')
            grid on
            xlabel('Displacement [mm]')
            ylabel('Damage [-]')
            ylim([0 1]);
        end

        function meshDamage(obj)
            obj.damageField.plot;
            title('Final damage distribution')
            shading interp
            colorbar
            clim([0 1])
        end

        function plotQR(obj)
            figure()
            plot(obj.rMax,obj.qMax)
            title('Hardening parameter (q-r)')
            xlabel('Hardening function (q) [$\sqrt{GPa}$]','Interpreter','latex')
            ylabel('Internal variable (r) [$\sqrt{GPa}$]','Interpreter','latex')
            hold off
        end
        
        function plotEnergy(obj)
            figure()
            plot(obj.displacement,obj.energy)
            title('Energy - displacement')
            xlabel('Displacement [mm]')
            ylabel('Energy [J]')
        end

        function plotIterations(obj)
            figure()
            plot(obj.iter)
            title('Iterations needed')
            xlabel('Step [-]')
            ylabel('Iterations')
        end

    end
    
end