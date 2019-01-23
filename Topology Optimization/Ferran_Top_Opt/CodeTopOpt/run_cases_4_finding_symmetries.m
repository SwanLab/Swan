function run_cases_4_finding_symmetries
cluster = 0;
opencluster(cluster)

Vfrac = 0.6;


nlevel = 7;
firstmesh = 1;
pint = 0; % 0 == no pintar, 1 == pintar
penalty = 1;
lambda0 = 0;


%[~,phi_all,theta_all,numsnapshots] = vademecum_spehere_mesh(nlevel,pint);

n_cuadrantes = 4;
phi_desv = pi/12;
theta_desv = pi/12;

phi_cuadr = [-pi/4 -3*pi/4 3*pi/4 pi/4]';
phi_minus = -phi_desv*ones(4,1);
phi_plus = phi_desv*ones(4,1);

phi_all = [phi_cuadr+phi_plus,phi_cuadr+phi_minus]';
phi_all = phi_all(:);
theta_all = pi/2-theta_desv*ones(size(phi_all));
strain = phitheta2strain(phi_all,theta_all);
A = zeros(3,3,size(strain,1));
for i = 1:size(strain,1)
A(:,:,i) = strain(i,:)'*strain(i,:);
end

nsnapshots = length(phi_all);



%nsnapshots = numsnapshots(end);

for isnapshot = 8:nsnapshots
isnapshot
%make_regula_falsi(isnapshot,phi_all,theta_all,Vfrac,firstmesh,penalty)
make_augmented(isnapshot,phi_all,theta_all,Vfrac,firstmesh,penalty,lambda0);

end



end


function strain = phitheta2strain(phi,theta)

stra_theta = theta;
stra_phi = phi;
stra_x = sin(stra_theta).*cos(stra_phi);
stra_y = sin(stra_theta).*sin(stra_phi);
stra_xy = cos(stra_theta);
strain = [stra_x stra_y stra_xy];

end


function opencluster(cluster)

if cluster == 1
addpath(genpath('/home/aferrer/matlab-hpc-interactive'))
addpath(genpath('/home/aferrer/MicrostructureTopologicalDerivative'))
nLAB = 412;
openmatlab(nLAB)
end

if cluster == 0
addpath(genpath('/home/aferrer/Documents/Doctorat/Tesi/MicroEscalaDT'))
end



end


function [VOL,cost,matCh] = make_augmented(isnapshot,phi_all,theta_all,Vfrac,firstmesh,penalty,lambda0)
    tipo_lambda = 'AUGMENTED';
    inicial = 1;
   
    [VOL,cost,matCh,problembsc,d_u,file_gid,iter,coordinates,element,dim_old,phifunct_n,g_nodal_n,gfunc_til_n,tstres,post,g_ortho,norm_g_n,vdisp,nbdata] = compute_Vol(isnapshot,penalty,phi_all,theta_all,Vfrac,firstmesh,lambda0,tipo_lambda,inicial);
    printinfo(1,problembsc,d_u,file_gid,iter,coordinates,element,dim_old,phifunct_n,g_nodal_n,gfunc_til_n,tstres,post,g_ortho,norm_g_n,vdisp,nbdata);
    
end




function make_regula_falsi(isnapshot,phi_all,theta_all,Vfrac,firstmesh,penalty)
    tipo_lambda = 'POTENCIAL';
    inicial = 1;
   
    lambda_min = 0;
    [VOL_max,cost_max,matCh_max] = compute_Vol(isnapshot,penalty,phi_all,theta_all,Vfrac,firstmesh,tipo_lambda,inicial);
    
    inicial = 0;
    lambda_max =  27;

    [VOL_min,cost_min,matCh_min] = compute_Vol(isnapshot,penalty,phi_all,theta_all,Vfrac,firstmesh,tipo_lambda,inicial);
       
    
    res = [VOL_max - Vfrac, VOL_min - Vfrac];
    lambda = [lambda_min lambda_max];
    cost = [cost_max,cost_min];
%     C11 =  [matCh_max(1,1),matCh_min(1,1)];
%     C12 =  [matCh_max(1,2),matCh_min(1,2)];
%     C13 =  [matCh_max(1,3),matCh_min(1,3)];
%     C22 =  [matCh_max(2,2),matCh_min(2,2)];
%     C23 =  [matCh_max(2,3),matCh_min(2,3)];
%     C33 =  [matCh_max(3,3),matCh_min(3,3)];
%         
    iter = 2;
    while abs(res(iter)) > 1e-4
        %irun = cases_to_run(isnapshot);
        lambda_new = max(0,((VOL_max-Vfrac)*lambda_max - (VOL_min-Vfrac)*lambda_min)/(VOL_max - VOL_min));
        [VOL_new,cost_new,matCh_new] = compute_Vol(isnapshot,penalty,phi_all,theta_all,Vfrac,firstmesh,lambda_new,tipo_lambda,inicial);
        
        res(iter+1) = VOL_new - Vfrac;
        lambda(iter+1) = lambda_new;
        cost(iter+1) = cost_new;
%         C11(iter+1) =  matCh_new(1,1);
%         C12(iter+1) =  matCh_new(1,2);
%         C13(iter+1) =  matCh_new(1,3);
%         C22(iter+1) =  matCh_new(2,2);
%         C23(iter+1) =  matCh_new(2,3);
%         C33(iter+1) =  matCh_new(3,3);
        
%         figure(50)
%         plot(abs(res))
%         
%         
%         [~,orden] = sort(lambda);
%         figure(100)
%         plot(lambda(orden),res(orden)+Vfrac,'-+')
%         hold on
%         plot(lambda(end),res(end)+Vfrac,'+r')
%         hold off
%         
%         figure(101)
%         plot(res(orden)+Vfrac,C11(orden),'-+')
%         legend('C11')
%         figure(102)
%         plot(res(orden)+Vfrac,C12(orden),'-+')
%         legend('C12')
%         figure(103)
%         plot(res(orden)+Vfrac,C13(orden),'-+')
%         legend('C13')
%         figure(104)
%         plot(res(orden)+Vfrac,C22(orden),'-+')
%         legend('C22')
%         figure(105)
%         plot(res(orden)+Vfrac,C23(orden),'-+')
%         legend('C23')
%         figure(106)
%         plot(res(orden)+Vfrac,C33(orden),'-+')
%         legend('C33')
   
        
        
        if res(iter+1) >= 0
           
            VOL_max = VOL_new;
            lambda_min = lambda_new;           
        else
            VOL_min = VOL_new;
            lambda_max = lambda_new;           
        end
        abs(res)
        iter = iter + 1;
         
    end

[VOL_new,cost_new,matCh_new,problembsc,d_u,file_gid,iter,coordinates,element,dim_old,phifunct_n,g_nodal_n,gfunc_til_n,tstres,post,g_ortho,norm_g_n,vdisp,nbdata] = compute_Vol(isnapshot,penalty,phi_all,theta_all,Vfrac,firstmesh,lambda_new,tipo_lambda,inicial);
printinfo(1,problembsc,d_u,file_gid,iter,coordinates,element,dim_old,phifunct_n,g_nodal_n,gfunc_til_n,tstres,post,g_ortho,norm_g_n,vdisp,nbdata);

end

function  [VOL,cost_n,matCh,problembsc,d_u,file_gid,iter,coordinates,element,dim_old,phifunct_n,g_nodal_n,gfunc_til_n,tstres,post,g_ortho,norm_g_n,vdisp,nbdata] = compute_Vol(isnapshot,penalty,phi_all,theta_all,Vfrac,firstmesh,lambda,tipo_lambda,inicial)

        irun = isnapshot;
        file_path = ['/home/aferrer/Documents/Doctorat/Tesi/MicroEscalaDT/Cases4FindingSymmetries/',num2str(irun),'/'];
        %unix(['rm -r ',file_path]);
        if inicial 
        unix(['mkdir ',file_path]);
        end
        phi = phi_all(irun);
        theta = theta_all(irun);

        [VOL,cost_n,matCh,problembsc,d_u,file_gid,iter,coordinates,element,dim_old,phifunct_n,g_nodal_n,gfunc_til_n,tstres,post,g_ortho,norm_g_n,vdisp,nbdata] = MAIN_OPT3(Vfrac,lambda,penalty,phi,theta,tipo_lambda,'MIN_STIFF_INV',file_path,firstmesh);

end



function openmatlab(nLAB)

nLabImp_TRAINING = nLAB;

nLab = matlabpool('size');
if nLab>0
    disp(['Se utiliza la configuracion activa: ',get(findResource(),'configuration'),...
        ' de ',num2str(nLab),' laboratorios.']);
else
    %Recordar que para pocos elementos no conviene paralelizar.
    if nLabImp_TRAINING==0
        matlabpool open
    elseif nLabImp_TRAINING>0
        matlabpool('open',nLabImp_TRAINING)
    end
end


end
