clear; clc; close all;

% --- Parameters ---
params.tiles_per_dim = [16; 16];          % 8x8 tiles per tiling
params.num_tilings = 4;                   % number of tilings
params.n_actions = 3;                     % {-1, 0, 1}
params.pos_range = [-1.2, 0.7];
params.vel_range = [-0.07, 0.07];
params.state_range = [params.pos_range; params.vel_range];  % 2 x 2
params.initial_state = [-0.5, 0];
params.n_features = prod(params.tiles_per_dim) * params.num_tilings * params.n_actions;

% SARSA(λ) hyperparameters
params.alpha = 0.2 / params.num_tilings;
params.gamma = 1.0;
params.lambda = 0.95;
params.epsilon = 0.1;
params.n_episodes = 4000;

% Actor-critic hyperparameters
params.alpha_c = 0.05/params.num_tilings;    % critic learning rate
params.alpha_a = 0.05/params.num_tilings;    % actor learning rate

% Epsilon decay strategy
params.epsilon_strategy = 'decay';   % or 'linear_decay'
params.epsilon_decay = 0.995;
%params.min_epsilon = 0.001;

% Plot resolution and gif frequency
params.pixels1 = 300;
params.pixels2 = 300;
params.gifFreq = 100;
params.max_steps = 4000;

% --- Offsets ---
params.offsets = generateDiagonalOffsets(2, params.num_tilings);  % 2D case

% --- Environment ---
env = MountainCarEnv(params);

% --- Old Versions ---
%agent = Agent(env, policyFunction, @getActiveTiles, params);
%[w_c, weights] = agent.ActorCritic();
%weights = agent.SARSA();

% --- Active Tiles ---
activeTiles = ActiveTiles(params);

% --- Policy ---
policyType = 'epsilonGreedy';
%policyType = 'softmax';
policyFunction = createPolicyFunction(policyType);

% --- Agent ---
agent = Sarsa(params.n_features);
%agent = Q_learning(params.n_features, activeTiles);
%agent = ActorCritic(params);

% --- Training ---
weights = TD_lambda(env, policyFunction, activeTiles, params, agent);
%[~, weights] = ActorCriticFunc(env, policyFunction, activeTiles, params, agent);
%% --- Visualization ---

% Resolution of plots
Pixels1 = params.pixels1;  % Position resolution
Pixels2 = params.pixels2;  % Velocity resolution

% Plot learned Q-value surface
plotQSurface2D(weights, params, Pixels1, Pixels2, activeTiles); 

% Plot learned policy
plotPolicyColorMap(weights, params, activeTiles);

maxSteps = params.max_steps;

% Plot single episode trajectory (greedy policy)
initialState = [-0.5; 0];
plotEpisodeTrajectory2D(weights, initialState, maxSteps, params, env, activeTiles)

% Plot multiple episode trajectories
nTrajectories = 10;
initialState = zeros(nTrajectories, 2);
for i = 1:nTrajectories
    pos = -0.8 + (0.2 + 0.8) * rand();
    vel = -0.05 + (0.02 + 0.05) * rand();
    initialState(i, :) = [pos, vel];
end
%initialState = [-0.5, 0];  % Starting at default initial state
plotMultipleTrajectories2D(weights, initialState, maxSteps, params, env, activeTiles);
