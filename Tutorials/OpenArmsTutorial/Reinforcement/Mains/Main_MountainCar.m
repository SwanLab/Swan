clear; clc; close all;

% --- Parameters ---
params.tiles_per_dim = [16; 16];          % 8x8 tiles per tiling
params.num_tilings = 4;                 % number of tilings
params.n_actions = 3;                   % {-1, 0, 1}
params.pos_range = [-1.2, 0.7];
params.vel_range = [-0.07, 0.07];
params.state_range = [params.pos_range; params.vel_range];  % 2 x 2
params.n_features = prod(params.tiles_per_dim) * params.num_tilings * params.n_actions;

% SARSA(λ) hyperparameters
params.alpha = 0.2 / params.num_tilings;
params.gamma = 1.0;
params.lambda = 0.9;
params.epsilon = 0.1;
params.n_episodes = 4000;

% Epsilon decay strategy
params.epsilon_strategy = 'decay';   % or 'linear_decay'
params.epsilon_decay = 0.995;
%params.min_epsilon = 0.001;

% Plot resolution and gif frequency
params.pixels1 = 300;
params.pixels2 = 300;
params.gifFreq = 100;

% --- Offsets ---
params.offsets = generateDiagonalOffsets(2, params.num_tilings);  % 2D case

% --- Environment ---
env = createMountainCarEnv(params);

% --- Train using SARSA(λ) ---
tic
weights = sarsa_lambda(env, @policyFunction, @getActiveTiles, params);
toc
save('w_reference.mat', 'weights');
%% --- Visualization ---

% Resolution of plots
Pixels1 = params.pixels1;  % Position resolution
Pixels2 = params.pixels2;  % Velocity resolution

% Plot learned Q-value surface
plotQSurface2D(weights, params, Pixels1, Pixels2); 

% Plot learned policy
plotPolicyColorMap(weights, params, @getActiveTiles);

% Plot single episode trajectory (greedy policy)
initialState = [-0.5; 0];  % Starting at default initial state
maxSteps = 1000;
plotEpisodeTrajectory2D(weights, initialState, maxSteps, params, env, @getActiveTiles);

function [a, epsilon] = policyFunction(s, w, epsilon, params)
    [a, epsilon] = epsilonGreedy(s, w, epsilon, params, @getActiveTiles);
end

