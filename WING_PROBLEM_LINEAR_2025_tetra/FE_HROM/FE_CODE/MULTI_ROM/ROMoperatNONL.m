function HROMVAR = ROMoperatNONL(DATAROM,BASES,DATA_REFMESH,DATAIN,HROMVAR) 

if nargin == 0
    load('tmp.mat')
end

% 1.  Reduced-order B matrix (all integration points)
% \BdomRED{e} =  \Bdom{e}  \BasisUdef{e}{}
load(DATAIN.NAME_WS_MODES,'Bdom')
HROMVAR.BdomRED = Bdom*DATAROM.BasisUdef; 

% 2. \BdomREDg{e} = \selectDOM{e} \BdomRED{e},
%  Reduced-order B matrix (selected integration points)
indPoints = small2large(HROMVAR.setPoints,HROMVAR.nstrain) ;
HROMVAR.BdomREDg = HROMVAR.BdomRED(indPoints,:) ; 

% 3. Reaction-def. disp. matrix 
% \Hqr{e}{} = \BasisUdef{e^T}{f}\BasisRdef{e}{f},
HROMVAR.Hqr = DATAROM.Hqr ;  % See BeamStiffnessMatrix.m 

% 4. Reaction-rig.b. disp, matrix 
%  \HqrRB{e}{} = \BasisUdef{e^T}{f}\BasisRrb{e}{f}
HROMVAR.HqrRB = DATAROM.HqrRB ;  % See BeamStiffnessMatrix.m 

% 5. Geometric matrix 
% \gRB{e}{} = \BasisUrb{e^T}{f}  \BasisRrb{e}{f},  % See ExtForcesROM_operators
HROMVAR.gRB= DATAROM.gRB ;


% 6. Compatibility matrix (def.)
% \Tcomp{e}{} =   \BasisRdef{e^T}{f} \DrotINTFini{e}{}   \BasisINTF{e}{},  %%%  See BeamStiffnessMatrix.m 
for i  =1:length(DATAROM.Tcomp)
    HROMVAR.Tcomp{i} = DATAROM.Tcomp{i}' ;
end
HROMVAR.Tcomp = cell2mat(HROMVAR.Tcomp ) ; 

% 6. Compatibility matrix (.)
% \TcompRB{e}{} =   \BasisRdef{e^T}{f} \DrotINTFini{e}{}   \BasisINTF{e}{},  %%%  See BeamStiffnessMatrix.m 
HROMVAR.TcompRB{1} = DATAROM.TcompRB{1} ; 
HROMVAR.TcompRB{2} = DATAROM.TcompRB{2} ; 
HROMVAR.TcompRB = cell2mat(HROMVAR.TcompRB ) ; 
%7. Coarse-scale B-matrix 
% \BCdom{e}{} =  \BdomREDg{e}  \Hqr{e^{-T}}\Tcomp{e}{} 
HROMVAR.BCdom =  HROMVAR.BdomREDg*(HROMVAR.Hqr'\HROMVAR.Tcomp) ;  

 
% REDUCED-ORDER OPERATORS 
% 
% 
%   \begin{equation*}
%  \label{eq.53**3q3}
%  \begin{split}
%   \hspace{0.5cm} \hspace{0.5cm}
%   \gRB{e}{} = \BasisUrb{e^T}{f}  \BasisRrb{e}{f}, \hspace{0.5cm} 
%   \Tcomp{e}{} =   \BasisRdef{e^T}{f} \DrotINTFini{e}{}   \BasisINTF{e}{}, \hspace{0.5cm}
%   \TcompRB{e}{} =    \BasisRrb{e^T}{f} \DrotINTFini{e}{}   \BasisINTF{e}{}, \hspace{0.5cm} 
%   \BCdom{e}{} =  \BdomREDg{e}  \Hqr{e^{-T}}\Tcomp{e}{} 
%   \end{split}
%  \end{equation*} 

% 



