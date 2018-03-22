        function checkDerivative(obj)
            obj.preProcess;
            Msmooth = obj.filter.Msmooth;
            x0 = obj.x;
            % Initialize function
            epsi = 1e-6;
            %initial
            compliance0 = ShFunc_Compliance(obj.settings);
            % compliance0.h_C_0 = 1;
            volume0 = ShFunc_Volume(obj.settings);
            perimeter0 = ShFunc_Perimeter(obj.settings);
            %new
            compliance = ShFunc_Compliance(obj.settings);
            % compliance.h_C_0 = 1;
            volume = ShFunc_Volume(obj.settings);
            perimeter = ShFunc_Perimeter(obj.settings);
            
            obj.incremental_scheme.update_target_parameters(1,compliance0,volume0,perimeter0);
            obj.incremental_scheme.update_target_parameters(1,compliance,volume,perimeter);
            %evaluate initial
            compliance0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            volume0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            perimeter0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            
            compliance0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            volume0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            perimeter0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            
            compliance.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            volume.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            perimeter.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            
            nnod = length(compliance0.gradient);
            g = zeros(nnod,1);
            gp = zeros(nnod,1);
            gv = zeros(nnod,1);
            for inode = 1:nnod
                if mod(inode,100) == 0
                    disp(strcat('Node: ',int2str(inode)));
                end
                xnew = x0;
                xnew(inode) = xnew(inode)-epsi;
                compliance.computef(xnew,obj.topOpt_params,obj.interpolation,obj.filter);
                volume.computef(xnew,obj.topOpt_params,obj.interpolation,obj.filter);
                perimeter.computef(xnew,obj.topOpt_params,obj.interpolation,obj.filter);
                g(inode) = (compliance0.value-compliance.value)/epsi;
                gv(inode) = (volume0.value-volume.value)/epsi;
                gp(inode) = (perimeter0.value-perimeter.value)/epsi;
            end
            fprintf('Relative error Volume: %g\n',obj.error_norm_field(gv,volume0.gradient,Msmooth));
            fprintf('Relative error Perimeter: %g\n',obj.error_norm_field(gp,perimeter0.gradient,Msmooth));
            fprintf('Relative error Compliance: %g\n',obj.error_norm_field(g,compliance0.gradient,Msmooth));
        end
        function enorm = error_norm_field(obj,gp,gp0,Msmooth)
            
            enodal = gp0 - gp;
            enorm = (enodal'*Msmooth*enodal)/(gp'*Msmooth*gp);
            
        end