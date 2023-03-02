% BPOPT Create a New Problem Instance
classdef bp_create < handle
    properties (Access = public)
        x
    end
    properties (Access = private)
    end

    methods (Access = public)
        function obj = bp_create(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.loadParameters();
            obj.selectProblem();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
        end

        function loadParameters(obj)
        % initial merit function weighting of constraint residuals
        obj.x.nu = 10;

        % initial barrier term
        obj.x.mu = 0.1;

        % update mu? - used for testing
        obj.x.mu_update = true;

        % line search criteria
        %  1 = reduction in merit function
        %  2 = simple clipping
        %  3 = filter method
        obj.x.line_search = 2;

        % Solving for z
        % 1 = update from direct solve approach
        % 2 = update explicitly from z(i) = mu / x(i)
        obj.x.z_update = 1;

        % reduced or full matrix inversion
        %  1 = condensed, symmetric matrix
        %  2 = full, unsymmetric matrix
        obj.x.matrix = 1;

        % show contour plots
        obj.x.contour = true;

        % maximum iterations
        obj.x.maxiter = 100;

        % Debug printing level (0-10)
        obj.x.idebug = 0;

        % Initialize slack variables to equation residuals (true)
        %   or to near zero 0.01 (false)
        obj.x.slack_init = true;

        % update tau
        obj.x.tau_max = 0.01;

        % solution error tolerance
        obj.x.e_tol = 1e-6;
        obj.x.k_e = 10;
        obj.x.k_soc = 0.99;

        % check for new barrier problem
        obj.x.k_mu  = 0.2; % (0,1)

        % not currently used
        obj.x.gamma_alpha = 0.05;
        obj.x.s_th = 1.1;  % > 1
        obj.x.s_phi = 2.3; % > 1
        obj.x.epsilon = 1; % > 0
        obj.x.eta_phi = 10^-4; % 0 < x < 0.5
        end

        function selectProblem(obj)
            if (isstr(name)),
                % Select custom problem
                % 0 = Custom model defined in APM
                obj.x.prob = 0;
                % application name
                obj.x.name = name;
                % server for obtaining problem information
                obj.x.server = 'http://xps.apmonitor.com';
                %x.server = 'http://localhost';
                % load application and get app name
                obj.x.app = apm_app(obj.x.server,name);
                % zero iterations - return stats only
                apm_option(obj.x.server,obj.x.app,'nlc.max_iter',0);
                % diaglevel must be at least 2
                apm_option(obj.x.server,obj.x.app,'nlc.diaglevel',4);    
                % line search with simple clipping
                obj.x.line_search = 2;
            else
                obj.x.prob = name;
                obj.x.name = 'test';
                obj.x.server = 'local';
                obj.x.app = 'test';
                
                switch (obj.x.prob)
                    case(1)
                        obj.x.nu = 1000;
                        obj.x.mu = 1e-1;
                        %x.z_update = 2;
                        %x.line_search = 2;
                    case(2)
                        obj.x.nu = 10;
                        obj.x.mu = 10;
                    case(3)
                        obj.x.nu = 10;
                        obj.x.mu = 1e-2;
                    case(4)
                        obj.x.nu = 10;
                        obj.x.mu = 10;
                        obj.x.maxiter = 50;
                        %x.mu_update = false;
                        %x.line_search = 2;
                        %x.z_update = 1;
                    case(5)
                        obj.x.nu = 10;
                        obj.x.mu = 10;
                    case(6)
                        obj.x.nu = 100;
                        obj.x.mu = 0.1;
                    case(7)
                        obj.x.nu = 100;
                        obj.x.mu = 0.001;
                    case(8)
                    case(9)
                    case(10)
                        %x.line_search = 2;
                        obj.x.slack_init = false;
                end
              end
        end
    end
end

%
% Inputs
%   name = application name to load
%
% Outputs
%   x = problem instance
function [x] = bp_create(name)

% initial merit function weighting of constraint residuals
x.nu = 10;

% initial barrier term
x.mu = 0.1;

% update mu? - used for testing
x.mu_update = true;

% line search criteria
%  1 = reduction in merit function
%  2 = simple clipping
%  3 = filter method
x.line_search = 2;

% Solving for z
% 1 = update from direct solve approach
% 2 = update explicitly from z(i) = mu / x(i)
x.z_update = 1;

% reduced or full matrix inversion
%  1 = condensed, symmetric matrix
%  2 = full, unsymmetric matrix
x.matrix = 1;

% show contour plots
x.contour = true;

% maximum iterations
x.maxiter = 100;

% Debug printing level (0-10)
x.idebug = 0;

% Initialize slack variables to equation residuals (true)
%   or to near zero 0.01 (false)
x.slack_init = true;

% update tau
x.tau_max = 0.01;

% solution error tolerance
x.e_tol = 1e-6;
x.k_e = 10;
x.k_soc = 0.99;

% check for new barrier problem
x.k_mu  = 0.2; % (0,1)

% not currently used
x.gamma_alpha = 0.05;
x.s_th = 1.1;  % > 1
x.s_phi = 2.3; % > 1
x.epsilon = 1; % > 0
x.eta_phi = 10^-4; % 0 < x < 0.5

if (isstr(name)),
    % Select custom problem
    % 0 = Custom model defined in APM
    x.prob = 0;
    % application name
    x.name = name;
    % server for obtaining problem information
    x.server = 'http://xps.apmonitor.com';
    %x.server = 'http://localhost';
    % load application and get app name
    x.app = apm_app(x.server,name);
    % zero iterations - return stats only
    apm_option(x.server,x.app,'nlc.max_iter',0);
    % diaglevel must be at least 2
    apm_option(x.server,x.app,'nlc.diaglevel',4);    
    % line search with simple clipping
    x.line_search = 2;
else
    x.prob = name;
    x.name = 'test';
    x.server = 'local';
    x.app = 'test';
    
    switch (x.prob)
        case(1)
            x.nu = 1000;
            x.mu = 1e-1;
            %x.z_update = 2;
            %x.line_search = 2;
        case(2)
            x.nu = 10;
            x.mu = 10;
        case(3)
            x.nu = 10;
            x.mu = 1e-2;
        case(4)
            x.nu = 10;
            x.mu = 10;
            x.maxiter = 50;
            %x.mu_update = false;
            %x.line_search = 2;
            %x.z_update = 1;
        case(5)
            x.nu = 10;
            x.mu = 10;
        case(6)
            x.nu = 100;
            x.mu = 0.1;
        case(7)
            x.nu = 100;
            x.mu = 0.001;
        case(8)
        case(9)
        case(10)
            %x.line_search = 2;
            x.slack_init = false;
    end
end