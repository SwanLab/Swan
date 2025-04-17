classdef PhaseFieldPlotter < handle
    
    properties (Access = private)
        damage   
        damageField
        displacement
        reaction
        energy
        iter
        costFun
    end
    
    methods (Access = public)
        
        function obj = PhaseFieldPlotter(cParams)
            obj.init(cParams)
            obj.meshDamage()
            obj.plotDamage()
            obj.plotForceDisplacement()
            obj.plotEnergies()
            obj.plotIterations()
            obj.plotCost()
        end
        
    end
        
    methods (Access = private)
        
        function init(obj,cParams)
            obj.damage = cParams.damage.maxValue;
            obj.damageField = cParams.damage.field;
            obj.displacement = cParams.displacement.value;
            obj.reaction = cParams.force;
            obj.energy = cParams.energy;
            obj.iter = cParams.iter;
            obj.costFun = cParams.cost;
        end
        
        function plotDamage(obj)
            figure()
            hold on
            plot(obj.displacement,obj.damage,'Color',"#0072BD")
            %obj.computeTheoreticalDamage()
            title('Damage-displacement diagram')
            grid on
            xlabel('Displacement [mm]')
            ylabel('Damage [-]')
            ylim([0 1]);
        end

        function computeTheoreticalDamage(obj)
            E = 210; nu = 0.3; l0 = 0.1; Gc = 5e-3;
            C22 = E/((1+nu)*(1-nu));
            e = obj.displacement;
            damageTheory = (C22.*e.^2)./((Gc/l0)+(C22.*e.^2));
            plot(obj.displacement,damageTheory);
        end
        
        function meshDamage(obj)
            obj.damageField.plot;
            title('Final damage distribution')
            shading interp
            colorbar
            clim([0 1])
        end
        
        function plotForceDisplacement(obj)
            figure()
            hold on
            plot(obj.displacement,obj.reaction)
            %obj.computeTheoreticalForce()
            title('Force-displacement diagram (Reaction force)')
            grid on
            xlabel('Displacement [mm]')
            ylabel('Force [kN]')
        end

        function computeTheoreticalForce(obj)
            E = 210; nu = 0.3; l0 = 0.1; Gc = 5e-3;
            C22 = E/((1+nu)*(1-nu));
            e = obj.displacement;
            damageTheory = (C22.*e.^2)./((Gc/l0)+(C22.*e.^2));
            C22phi = ((1-damageTheory).^2).*C22;
            sigma = C22phi.*e;
            plot(obj.displacement,sigma);
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

        function plotCost(obj)
            figure()
            hold on
            plot(obj.costFun(1,:))
            title('Cost Function')
            xlabel('Iteration [-]')
            ylabel('Energy [J]')
            hold off
        end

    end
    
end