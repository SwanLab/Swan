%% LPV model definition

% Model constants
a=-0.4;
c_l_beta_1=1.6790;
c_l_beta_2=2.3506;
c_m_beta_1=-0.3175;
c_m_beta_2=-0.3175;
b=0.135;
c_alfa=0.036;
m=12.387;
k_h=2844.4;
ro=1.225;
I_alfa=0.065;
c_h=27.43;
c_l_alfa=6.28;
c_m_alfa=(0.5+a)*c_l_alfa;
x_alfa=-0.3533-a;
c_l_beta=3.358;
c_m_beta=-0.635;

d=m*(I_alfa-m*x_alfa^2*b^2);
k_1=(I_alfa*k_h)/d;
k_2=(I_alfa*ro*b*c_l_alfa+m*x_alfa*b^3*ro*c_m_alfa)/d;
k_3=(-m*x_alfa*b*k_h)/d;
k_4=(-m*x_alfa*b^2*ro*c_l_alfa-m*ro*b^2*c_m_alfa)/d;
c_1=@(U)(I_alfa*(c_h+ro*U*b*c_l_alfa)+m*x_alfa*ro*U*b^3*c_m_alfa)/d;
c_2=@(U)(I_alfa*ro*U*b^2*c_l_alfa*(1/2-a)-m*x_alfa*b*c_alfa+m*x_alfa*ro*U*b^4*c_m_alfa*(1/2-a))/d;
c_3=@(U)(-m*x_alfa*b*c_h-m*x_alfa*ro*U*b^2*c_l_alfa-m*ro*U*b^2*c_m_alfa)/d;
c_4=@(U)(m*c_alfa-m*x_alfa*ro*U*b^3*c_l_alfa*(1/2-a)-m*ro*U*b^3*c_m_alfa*(1/2-a))/d;

g_3=(1/d)*(-I_alfa*ro*b*c_l_beta-m*x_alfa*b^3*ro*c_m_beta);
g_4=(1/d)*(m*x_alfa*b^2*ro*c_l_beta+m*ro*b^2*c_m_beta);

k_alfa_x2=@(alfa)2.82*(1-22.1*alfa+1315.5*alfa^2-8580*alfa^3+17289.7*alfa^4);

p_x2=@(alfa)(-m*x_alfa*b*k_alfa_x2(alfa))/d;
q_x2=@(alfa)m*k_alfa_x2(alfa)/d;

% LPV model
LPV = {...
    @(p)0     @(p)0     @(x)1     @(p)0     @(p)0; 
    @(p)0     @(p)0     @(x)0     @(p)1     @(p)0;
    @(p)-k_1  @(p)-k_2*p(1)^2-p_x2(p(2)) @(p)-c_1(p(1)) @(p)-c_2(p(1)) @(p)g_3*p(1)^2;
    @(p)-k_3  @(p)-k_4*p(1)^2-q_x2(p(2)) @(p)-c_3(p(1)) @(p)-c_4(p(1)) @(p)g_4*p(1)^2;
};

% parameter dependencies
dep = zeros([size(LPV) 2]);
dep(3,2,:) = [1 1];
dep(3,3,:) = [1 0];
dep(3,4,:) = [1 0];
dep(3,5,:) = [1 0];
dep(4,2,:) = [1 1];
dep(4,3,:) = [1 0];
dep(4,4,:) = [1 0];
dep(4,5,:) = [1 0];

n = 4;
