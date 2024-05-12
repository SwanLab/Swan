%**************************************************************************
% Shape Functional Associated to the Compliance with Volume Constraint 
%**************************************************************************
% DESCRIPTION
% Computes the shape functional
% 
% INPUT
% F:       pdetool load
% U:       pdetool solution
% volume:  current volume
% params:  topology optimization parameters struct
%
% OUTPUT
% sf: shape functional value
%
% HISTORY
% S. Amstutz     06/2009: code implementation.
% A.A. Novotny   06/2009: code updating.
% J-M.C. Farias  12/2010: code updating.
% D.E. Campeão   
% A.A. Novotny
% A. Romero  6/2018: code updating multimaterial.
%**************************************************************************

function [sf,Epot] = shfunc(F,U,volume,params,TOParams)

    penalty = params.penalty; max_vol = TOParams.max_vol;  auglag= params.auglag;
    volfrac = params.volfrac; voltarget = max_vol.*volfrac;
    sf = 0.5*dot(F,U) / TOParams.energy0 ; %normalizing compliance
    Epot = sf;
    if params.penalization == 1 % linear penalization
        sf = sf + sum( penalty .* volume(1:end-1) ./ max_vol ) ; 
    elseif params.penalization == 2 % exact quadratic penalization
        coef = volume(1:end-1) ./ voltarget;
        sf = sf + penalty .* ((coef <= 1)*0.0+(coef > 1).*(1.0-coef).^2);
    elseif params.penalization == 3 % augmented lagrangian multiplier
        coef = volume(1:end-1) ./ voltarget;
%         aux = max(coef-1.0, -penalty ./ auglag);
%         sf = sf + sum(penalty .* aux + (aux).^2 .* auglag / 2.0); 
         
        aux = penalty - auglag.*(1.0-coef); %testing equality constraints
        sf = sf + sum((1.0./(2.0.*auglag)) .* ((aux).^2 - penalty.^2)); %testing equality constraints
        
        
    end

end

