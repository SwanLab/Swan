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
            obj.displacement = cParams.displacement;
            obj.reaction = cParams.reaction;
            obj.energy = cParams.energy;
            obj.iter = cParams.iter;
            obj.costFun = cParams.cost;
        end
        
        function plotDamage(obj)
            figure()
            plot(obj.displacement,obj.damage,'Color',"#0072BD")
            title('Damage-displacement diagram')
            grid on
            xlabel('Displacement [mm]')
            ylabel('Damage [-]')
            ylim([0 1]);
            
            % quad = Quadrature.set(obj.mesh.type);
            % quad.computeQuadrature('QUADRATIC');
            % obj.materialPhaseField.computeMatIso(quad);
            % C = obj.materialPhaseField.material.C(2,2,1,5);
            % Gc = obj.materialPhaseField.Gc;
            % e = obj.displacementMat(step);
            % obj.damageTheory(step) = (C*e^2)/((Gc/obj.l0)+(C*e^2));
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

        function plotCost(obj)
            figure()
            hold on
            for n = 2:size(obj.costFun,2)
                if obj.costFun(2,n) == 0
                    plot([n-1, n],[obj.costFun(1,n-1), obj.costFun(1,n)],'b')
                elseif obj.costFun(2,n) == 1
                    plot([n-1, n],[obj.costFun(1,n-1), obj.costFun(1,n)],'r')
                elseif obj.costFun(2,n) == 2
                    plot([n-1, n],[obj.costFun(1,n-1), obj.costFun(1,n)],'k')
                end
            end
            title('Cost Function')
            xlabel('Iteration [-]')
            ylabel('Energy [J]')
            hold off
        end

    end
    
end