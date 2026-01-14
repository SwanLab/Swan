function   Bbar =   BBarMATRIXquad(B1,B2,Bbar1,Bbar2)
% BBar matrix
%dbstop('4')
if nargin == 0
    load('tmp.mat')
end

% Dilatational part  (standard)
% -----------------
B = [B1 0 
    0   B2 
    B2  B1
    0  0
    ];
 
Bdil = 1/3*[B1 B2
    B1 B2
    0   0
      B1  B2  ] ;
% Deviatoric part
Bdev = B- Bdil ;
%%% Bbar dilatational part
Bbar_dil = 1/3*[Bbar1 Bbar2
                Bbar1 Bbar2
                  0      0
                Bbar1 Bbar2               
                  ]  ;

            
% Bbar matrix 
Bbar = Bdev + Bbar_dil ; 





