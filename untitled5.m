eps   = sym("eps",[1 3]);
eps   = [1 0 0];
theta = sym("theta");

E  = 1;
nu = 0.3;
k  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E,nu,2);
mu = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E,nu);
l  = IsotropicElasticMaterial.computeLambdaFromShearAndBulk(mu,k,2);
C = zeros(3,3);
C(1,1)= 2*mu+l;
C(1,2)= l;
C(2,1)= l;
C(2,2)= 0;%2*mu+l;
C(3,3)= mu;
Reps = [(1+cos(2*theta))/2 , (1-cos(2*theta))/2 , (sin(2*theta))/2 ;
        (1-cos(2*theta))/2 , (1+cos(2*theta))/2 , (-sin(2*theta))/2  ;
        -sin(2*theta)       , sin(2*theta)      , (cos(2*theta))    ];

Rsig = [(1+cos(2*theta))/2 , (1-cos(2*theta))/2 , sin(2*theta)   ;
        (1-cos(2*theta))/2 , (1+cos(2*theta))/2 , -sin(2*theta)  ;
        (-sin(2*theta))/2  , (sin(2*theta))/2   , (cos(2*theta)) ];

energy = eps*Rsig*C*Reps*eps.';
dEnergy = diff(energy,theta);
solve(dEnergy==0,theta)