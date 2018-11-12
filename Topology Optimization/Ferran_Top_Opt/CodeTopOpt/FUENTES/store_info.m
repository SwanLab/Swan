function [vcost_theta] = store_info(i,cost,theta,kappa,lambda,h_C,constr,Energy,vol,vcost_theta)

    vcost_theta(i,1)=i;
    vcost_theta(i,2)=cost;
    vcost_theta(i,3)=theta;
    vcost_theta(i,4)=kappa;
    vcost_theta(i,5)=lambda;
    vcost_theta(i,6)=h_C;
    vcost_theta(i,7)=constr;
    vcost_theta(i,8)=Energy;
    vcost_theta(i,9)=vol;

end

