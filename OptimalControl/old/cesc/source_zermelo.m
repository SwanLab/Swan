function [states_dot] = source_zermelo(states, controls) %#codegen

% model interface created by falcon.m

% Extract states
x = states(1);
y = states(2);

% Extract controls
u = controls(1);

% Constants
V=30; % m/s

% ------------------------ %
% implement the model here %
% ------------------------ %

% implement state derivatives here
x_dot = V*cos(u) + 0.5*y;
y_dot = V*sin(u);
states_dot = [x_dot; y_dot];
% EoF
end