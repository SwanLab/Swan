function [kappa_opt] = kappa_frequency(dim,phifunct_n,g_nodal_n,norm_g,theta_n,...
    element,coordinatesn,coordinatesa,problembsc)
% cuenta el numero de cambios de signo de phi 

TOL=1e-4;
maxiter=100;
kappa_min = 0;
kappa_max = 1;
percent = 5; % en porcentaje, el numero con los que me quedo
npnod = dim.npnod;
kappas = zeros(npnod,1);

% lazo sobre todos los nodos de la malla
for i=1:npnod
    phi_a = phifunct_n(i);
    phi_b = g_nodal_n(i)/norm_g;
    xa = kappa_min; xb=kappa_max; 
    
    % estudia el cambio de signo de phi
    if (phi_a*phi_b < 0) 
        % metodo de biseccion
        iter=1;
        while ((xb-xa)>TOL && iter<=maxiter)
            x = (xb+xa)/2;
            %se deberia hacer nodal .... esta es muy cara
            [phi] = update_phifunc(theta_n,x,phifunct_n,g_nodal_n,norm_g);
            phi_i = phi(i);
            if (phi_a*phi_i < 0)
               phi_b = phi_i;
               xb = x; 
            elseif (phi_b*phi_i <= 0)
                phi_a = phi_i;
                xa = x;
            end
            iter=iter+1;
        end
        kappas(i) = x;
    else
        kappas(i) = kappa_min;
    end
end
vkappa=nonzeros(kappas);
nvkappa = size(vkappa,1);
if (nvkappa>0) 
    vkappa=sort(vkappa);
    imax = round(size(vkappa,1)*percent/100)+1;
    kappa_opt = vkappa(imax);
else
    kappa_opt = 1e-9;
end
vol = zeros(nvkappa,1);
for i=1:nvkappa
    kappa_i = vkappa(i);
    [phi_i] = update_phifunc(theta_n,kappa_i,phifunct_n,g_nodal_n,norm_g);
    [vol_omega] = cal_omega(dim.nelem,element,dim.ndime,dim.nnode,...
        coordinatesn,coordinatesa,problembsc.problemtype,phi_i,dim.npnod);
    vol(i)=vol_omega;
end

figure(1)
hold on;
x = linspace(0,1,nvkappa);
plot(vkappa,vol);
end

