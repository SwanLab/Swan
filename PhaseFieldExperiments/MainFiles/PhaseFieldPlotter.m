classdef PhaseFieldPlotter < handle
    
    properties (Access = private)
        damage   
        damageField
        displacement
        reaction
        energy
        iter
    end
    
    methods (Access = public)
        
        function obj = PhaseFieldPlotter(cParams)
            obj.init(cParams)
          %  obj.meshDamage()
            obj.plotDamage()
            obj.plotForceDisplacement()
            obj.plotEnergies()
            obj.plotIterations()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.damage = cParams.damage.maxValue;
         %   obj.damageField = cParams.damage.field;
            obj.displacement = cParams.displacement;
            obj.reaction = cParams.reaction;
            obj.energy = cParams.energy;
            obj.iter = cParams.iter;
        end
        
        function plotDamage(obj)
            figure()
            plot(obj.displacement,obj.damage,'Color',"#0072BD")
            title('Damage-displacement diagram')
            grid on
            xlabel('Displacement [mm]')
            ylabel('Damage [-]')
            ylim([0 1]);
        end
        
        function meshDamage(obj)
            figure()
            obj.damageField.plot;
            title('Final damage distribution')
            colorbar
            clim([0 1])
        end
        
        function plotForceDisplacement(obj)
            figure()
            plot(obj.displacement,obj.reaction)
            title('Force-displacement diagram (Reaction force)')
            grid on
            xlabel('Displacement [mm]')
            ylabel('Force [kN]')
        end
        
        function plotEnergies(obj)
            figure()
            hold on
            
            E = struct2cell(obj.energy);
            totE = zeros(size(E{1}));
            for i=1:length(E)
                totE = totE + E{i};
            end

            plot(totE,'Color',"#0072BD")
            plot(obj.energy.intE,'Color',"#D95319")
            plot(obj.energy.localDis,'Color',"#EDB120")
            plot(obj.energy.regDis,'Color',"#7E2F8E")
            plot(obj.energy.extWork)
            hold off
            title('Energy values at each step (Equation terms)')
            legend('Total energy', ...
                   'Internal Energy', ...
                   'Local surface energy', ...
                   'Non-local surface energy', ...
                   'External Work')
            xlabel('Step [-]')
            ylabel('Energy [J]')
        end

        function plotIterations(obj)
            figure()
            hold on
            plot(obj.iter.u,'xb')
            plot(obj.iter.phi,'r')
            plot(obj.iter.stag,'k')
            hold off
            title('Iterations needed')
            legend('U','phi')
            xlabel('Step [-]')
            ylabel('Iterations')
        end

    end
    
end