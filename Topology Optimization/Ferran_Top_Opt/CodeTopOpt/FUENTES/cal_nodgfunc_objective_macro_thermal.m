function [gfunc] =  cal_nodgfunc_objective_macro_thermal(igaus,stres,strain,phifunct,element,problembsc,dim,vol_void,d_u,fextpoint)
              
nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

kxx = ones(1,nelem).*element.material.k11;
kxy = ones(1,nelem).*element.material.k12;
kyy = ones(1,nelem).*element.material.k22;

[chi] = cal_caracteristic_function(igaus,phifunct,element,dim,problembsc,vol_void);
Ce = fourier_law_elas(chi,kxx,kxy,kyy);

work1 = zeros(1,nelem);
for istre=1:nstre
    work1 = work1 + stres(istre,:).*strain(istre,:);
end
        
detK = kxx.*kyy;
detKs = sqrt(detK);

alpha = 1./sqrt(element.material.k11);
beta = 1./sqrt(element.material.k22);

gamma  = element.material.opt_epsi;

P_menos = polarization_tensor(gamma,alpha,beta,nstre,nelem);
P_mas = polarization_tensor(1/gamma,alpha,beta,nstre,nelem);

gamma = 1;
delta_menos = gamma;
delta_mas = 1/gamma;

g_mas = zeros(1,nelem);
g_menos = zeros(1,nelem);
for istre = 1:nstre
    for kstre = 1:nstre
        for jstre = 1:nstre
            g_mas = g_mas + strain(istre,:).*shiftdim(Ce(istre,kstre,:).*P_mas(kstre,jstre,:),1).*strain(jstre,:);
            g_menos = g_menos + strain(istre,:).*shiftdim(Ce(istre,kstre,:).*P_menos(kstre,jstre,:),1).*strain(jstre,:);
            
        end
    end
end

upg = interpol(d_u,element,dim,problembsc);
fpg = interpol(fextpoint,element,dim,problembsc);


g_menos  = -detKs.*g_menos + (1 - delta_menos)*fpg.*upg;
g_mas = -detKs.*g_mas + (1 - delta_mas)*fpg.*upg;


% epsi  = element.material.opt_epsi;
% mu = element.material.mu;
% lambda = element.material.lambda;
% coeff = (lambda+3*mu)/(lambda+mu);
% 
% c1 = (epsi-1)/(coeff*epsi+1)*(coeff+1)/2;
% c2 = (epsi-1)*(coeff-2)/(coeff + 2*epsi -1);
% alpha1 = 2*c1;
% alpha2 = c1*c2;
% 
% d1 = (1-epsi)/(coeff+epsi)*(coeff+1)/2;
% d2 = (1-epsi)*(coeff-2)/(coeff*epsi + 2 -epsi);
% beta1 = 2*d1;
% beta2 = d1*d2;
% 
% g_menos = alpha1*work1 + alpha2*work2  ;
% g_mas = beta1*work1 + beta2*work2 ;

gfunc(igaus,:) =   vol_void.*(g_mas) + (1-vol_void).*(g_menos);



% gfunc(igaus,:)= -g_menos;
% outside= (phigp > 0);
% if any(outside)
%     gfunc(igaus,outside)= g_mas(igaus,outside);
% end

 

end 


function P = polarization_tensor(gamma,alpha,beta,nstre,nelem)
P = zeros(nstre,nstre,nelem);
pfactor = 0.5*(1-gamma).*alpha.*beta;
P(1,1,:) = pfactor.*(alpha + beta)./(alpha + gamma.*beta);
P(2,2,:) = pfactor.*(alpha + beta)./(beta + gamma.*alpha);
end