% TP model transformation based controller design for the parallel-type double inverted pendulum
% Szabolcs Nagy, Zoltan Petres, Peter Baranyi
% FUZZ-IEEE 2008 p1374-1380

m_k = 1.0; % kg
m_1 = 0.3; % kg
m_2 = 0.1; % kg
g   = 9.8; % m/s^2
l_1 = 0.6; % m
l_2 = 0.2; % m

F = @(p)(1-3/4*cos(p(1))^2)*m_1;
G = @(p)(1-3/4*cos(p(2))^2)*m_2;
H = @(p)4/3*(m_k+F(p)+G(p));

% state vector: pos, a, b, dpos, da, db
% parameter vector: a, b, da, db
LPV = {...
  @(p)0 @(p)0 @(p)0 @(p)1 @(p)0 @(p)0 @(p)0;
  @(p)0 @(p)0 @(p)0 @(p)0 @(p)1 @(p)0 @(p)0;
  @(p)0 @(p)0 @(p)0 @(p)0 @(p)0 @(p)1 @(p)0;
  @(p)0 @(p)-m_1*g*sinc(p(1)/pi)*cos(p(1))/H(p)                  @(p)-m_2*g*sinc(p(2)/pi)*cos(p(2))/H(p)                  @(p)0 @(p)4/3*m_1*l_1*p(3)*sin(p(1))/H(p)            @(p)4/3*m_2*l_2*p(4)*sin(p(2))/H(p)            @(p)4/3/H(p);
  @(p)0 @(p)(m_k+m_1+G(p))*g*sinc(p(1)/pi)/l_1/H(p)              @(p)3/4*m_2*g*sinc(p(2)/pi)*cos(p(2))*cos(p(1))/l_1/H(p) @(p)0 @(p)-m_1*p(3)*sin(p(1))*cos(p(1))/H(p)         @(p)-m_2*l_2*p(4)*sin(p(2))*cos(p(1))/l_1/H(p) @(p)-cos(p(1))/l_1/H(p);
  @(p)0 @(p)3/4*m_1*g*sinc(p(1)/pi)*cos(p(1))*cos(p(2))/l_2/H(p) @(p)(m_k+m_2+F(p))*g*sinc(p(2)/pi)/l_2/H(p)              @(p)0 @(p)-m_1*l_1*p(3)*sin(p(1))*cos(p(2))/l_2/H(p) @(p)-m_2*p(4)*sin(p(2))*cos(p(2))/H(p)         @(p)-cos(p(2))/l_2/H(p);
};

%Parameter relation tensor
dep = zeros([size(LPV) 4]);
dep(4,2,:) = [1 1 0 0];
dep(4,3,:) = [1 1 0 0];
dep(4,5,:) = [1 1 1 0];
dep(4,6,:) = [1 1 0 1];
dep(4,7,:) = [1 1 0 0];
dep(5,2,:) = [1 1 0 0];
dep(5,3,:) = [1 1 0 0];
dep(5,5,:) = [1 1 1 0];
dep(5,6,:) = [1 1 0 1];
dep(5,7,:) = [1 1 0 0];
dep(6,2,:) = [1 1 0 0];
dep(6,3,:) = [1 1 0 0];
dep(6,5,:) = [1 1 1 0];
dep(6,6,:) = [1 1 0 1];
dep(6,7,:) = [1 1 0 0];

n = 6;
