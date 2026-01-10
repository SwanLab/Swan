clc
clear all
% INPUTS
% -------------------
format long g
NameFileMesh = 'malla1.msh';   % Name of the mesh
DATA.NAME_FILE_STORE = 'DATA_MASSK' ;  % Rootname
nameWS = [DATA.NAME_FILE_STORE,'_',NameFileMesh,'.mat'] ;
nMODES = 6 ; % Number of modes to be included in the modal representation
DATA.nMODESplot =nMODES ;
DATA.nsteps = 500 ;
DATA.ntimesPERIOD = 20 ;
DATA.imode_EVOLUTION =0 ; % Only plot one mode per simulation
DATA.xiDmin = 0.02; % Minimum damping ratio
% END INPUTs
% Input data generated with mainELASTOSTATIC.m
load(nameWS,'K','M','d','DOFl','DOFr','COOR','CN','TypeElement','posgp','NameFileMesh')


%%% Mass and stifness matrices (unrestricted DOFs)
Kall = K ; Mall = M ;
K = Kall(DOFl,DOFl) ;
M = Mall(DOFl,DOFl) ;

%%%% Determination of natural frequencies and modes
% -------------------------------------------------
[BasisU FREQ] = UndampedFREQ(M,K,nMODES)  ;

%% Plotting natural modes
GidPostProcessModes(COOR,CN,TypeElement,BasisU,posgp,NameFileMesh,DATA,DOFl);

%% Initial conditions
% ---------------------
dINI = d(DOFl) ;
vINI =zeros(size(dINI));
% Therefore
qDini = zeros(nMODES,1) ;
qINI = BasisU'*dINI ;
figure(1)
hold on
xlabel('Mode')
ylabel('(Amplitude)')
bar(abs(qINI))

%%% Evolution in time
LargestPeriod = 2*pi/FREQ(1) ;
T = DATA.ntimesPERIOD*LargestPeriod ;
t = linspace(0,T,DATA.nsteps) ;
DISP = zeros(length(d),DATA.nsteps) ;

%%%% Damping ratio
xiD = CalcDampingRatios(DATA,FREQ) ;
% Damped frequencies 
FREQd = FREQ.*sqrt(1-xiD.^2) ; 


%%%%%

if DATA.imode_EVOLUTION == 0
    for imode=1:size(BasisU,2)
        
        % Exponential term 
        termEXP = exp(-xiD(imode)*FREQ(imode)*t) ; 
        % Cosine term 
        termCOS = qINI(imode)*cos(FREQd(imode)*t) ; 
         % Sine term 
         factor = (qDini(imode) + xiD(imode)*FREQ(imode)*qINI(imode))/FREQd(imode) ; 
        termSIN =  factor*sin(FREQd(imode)*t);         
        qi = termEXP.*(termCOS + factor*termSIN)  ;
        qi = repmat(qi,length(DOFl),1) ;
        qi = bsxfun(@times,qi,BasisU(:,imode)) ;
        DISP(DOFl,:) = DISP(DOFl,:) + qi ;
    end
else
    error('Option not implemented')
%     imode = DATA.imode_EVOLUTION ;
%     qi = qINI(imode)*cos(FREQ(imode)*t) ;
%     qi = repmat(qi,length(DOFl),1) ;
%     qi = bsxfun(@times,qi,BasisU(:,imode)) ;
%     DISP(DOFl,:) = DISP(DOFl,:) + qi ;
end

GidPostProcessDynamic(COOR,CN,TypeElement,DISP,'DYNAMICdamp',posgp,NameFileMesh,t);