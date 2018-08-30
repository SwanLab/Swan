function [kappa_opt,liter,vcostl,remesh] = kappa_best(phifunct_l,g_nodal_n,theta_n,theta_t,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,alpha,h_C_0,Perim0,lambda_n,penalty_n,eta_n,iter,costfunc_n,vol_n,file_name,kappa_ini,Msmooth,Stiff_smooth,emass,flag_change_micro_ini,Group,file_gid,file_write,u0)
                                                     
algorithm_update = problembsc.algorithm_update;
lmax = problembsc.nkappa;
kend = element.material.kappa_min;
vcostl = zeros(lmax,1);
remesh = 0;


switch problembsc.algorithm_update
    case 'AMSTUTZ'
        kappa = linspace(1,kend,lmax);
    case 'BCN'
        kappa = linspace(cos(theta_n),kend,lmax);
    case 'BCN_gphi_SI'    
        kappa = linspace(1,kend,lmax);
    case 'BCN_gphi_IM'
        kappa = linspace(1/(norm_g_n*sin(theta_n)),kend,lmax);
    case 'BCN_phig_SI'
%         int1_mod = (1-sin(theta_n)) - kend;
%         int2_mod = (2-(1+sin(theta_n)));
%         int_t_mod = int1_mod + int2_mod;
%         prop1 = int1_mod/(int_t_mod);
%         prop2 = int2_mod/(int_t_mod);
%         int1 = linspace(1-sin(theta_n),kend,round(prop1*lmax));
%         int2 = linspace(2,1+sin(theta_n),round(prop2*lmax));
%         kappa = [int1 int2];
       
        kappa = linspace(1-sin(theta_n),kend,lmax);
        
    case 'BCN_phig_EX'
        kappa = linspace((1-sin(theta_n))/sin(theta_n),kend,lmax);
        
end


%kappa = linspace(kini,1,lmax); 
cost_seg = zeros(length(kappa),1);
theta_seg = zeros(length(kappa),1);
vol_seg = zeros(length(kappa),1);
flag_up_date_lambda = 0;
for ikappa = 1:length(kappa)
    kappa_l = kappa(ikappa);
    [phifunct_l1] = update_phifunc(theta_n,kappa_l,phifunct_l,g_nodal_n,norm_g_n,remesh,phifunct_l,problembsc.algorithm_update);
    [costk,thetak,~,volk] = equilibrium_update(phifunct_l1,alpha,h_C_0,Perim0,lambda_n,penalty_n,eta_n,theta_n,theta_t,element,problembsc,fixnodes,coordinates,fext,dim,Msmooth,Stiff_smooth,emass,flag_up_date_lambda,flag_change_micro_ini,Group,iter,file_gid,file_write,u0);
    %[costk,thetak,volk] = call_new_cost(length(kappa)-ikappa,kappa(ikappa),theta_n,phifunct_l,g_nodal_n,norm_g_n,element,dim,problembsc,fixnodes,coordinates,fext,emass);
    cost_seg(ikappa,1) = costk;
    theta_seg(ikappa,1) = thetak;
    vol_seg(ikappa,1) = volk;
end

[~,i_min] = min(cost_seg);
kappa_opt = kappa(i_min);
liter = 1;

% figure(23)
% 
% [~,hLine1,hLine2] = plotyy([kappa',linspace(0,1,length(kappa))'],[cost_seg,costfunc_n*ones(length(kappa)',1)],[kappa',linspace(0,1,length(kappa))'],[theta_seg,theta_n*ones(length(kappa)',1)])
% set(hLine1(1),'Marker','+')
% set(hLine2(1),'Marker','+')
% set(hLine1,'Color','b')
% set(hLine2,'Color','k')
% 

drawnow
figure(32)
plot(kappa',cost_seg,'+-')
 

end


