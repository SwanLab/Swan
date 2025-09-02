classdef SuperformulaFunctionality
    
    methods (Access = public)

        function rad = calculate(obj, phi, a, b, m, n1, n2, n3)

            rad = (abs(cos(m.*phi./4)./a).^n2 + abs(sin(m.*phi./4)./b).^n3).^(-1/n1);
        
        end

        function rad = find_min_radius(obj, a, b, m, n1, n2, n3)

            rad = min(obj.calculate(linspace(0, 2*pi, 1000), a, b, m, n1, n2, n3));
        
        end
        
        function rad = find_max_radius(obj, a, b, m, n1, n2, n3)
        
            rad = max(obj.calculate(linspace(0, 2*pi, 1000), a, b, m, n1, n2, n3));
        
        end

        function area = compute_area(obj, a, b, m, n1, n2, n3)
        
            theta_vec = linspace(0, 2*pi, 1e3);
            rad_vec = obj.calculate(theta_vec, a, b, m, n1, n2, n3);
            theta_diff = theta_vec(2) - theta_vec(1);
            area = 0.5 * sum(rad_vec.^2 .* theta_diff);
        
        end

        function is_valid = evaluate(obj, gPar, r_min, r_max)
        
            is_valid = true;
            tol = 1e-14;
        
            a = gPar.semiHorizontalAxis;
            b = gPar.semiVerticalAxis;   
            
            m = gPar.m;
            n1 = gPar.n1;
            n2 = gPar.n2;
            n3 = gPar.n3;
            
            periodicity_residue = abs(abs(cos(pi*m/2))^n2+ a ^n2/(b^n3)*abs(sin(pi*m/2))^n3 - 1);
            large_size_residue = max(obj.find_max_radius(a, b, m, n1, n2, n3) - r_max, 0);
            small_size_residue = max(r_min - obj.find_min_radius(a, b, m, n1, n2, n3), 0);
        
            if periodicity_residue > tol || large_size_residue > 0 || small_size_residue > 0
                is_valid = false;
            end
        
        end

        function a = find_closing_a(obj, gPar, r_target)
        
            b = gPar.semiVerticalAxis;
            m = gPar.m;
            n1 = gPar.n1;
            n2 = gPar.n2;
            n3 = gPar.n3;
            obj_fun = @(a) obj.find_max_radius(a, b, m, n1, n2, n3) - r_target;

            try
                a = fzero(obj_fun, [1e-16, 10]);
                if a < 1e-6 || obj.find_min_radius(a, b, m, n1, n2, n3) < 1e-2
                    a = -1;
                end
            catch
                a = -1;
            end
        
        end

    end
end