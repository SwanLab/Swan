function VAR = CheckInSignDamageModel(VAR,DATA) 
if nargin == 0
    load('tmp.mat')
end

if isfield(DATA.DATA_evaluateTAU_and_DER,'signqNONrel_training') 
    qLIN = VAR.DISP(1) ; 
    qNON = VAR.DISP(2) ; 
    if sign(qLIN*qNON) ~= DATA.DATA_evaluateTAU_and_DER.signqNONrel_training
        fprintf('Flipping sign of qNON\n');
        qNON = -qNON ; 
    end
    VAR.DISP(2) = qNON ;  
end



