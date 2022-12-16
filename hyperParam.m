%% predefined hyper parameters for the algorithm
% below the used units are provided. Using known already units.
% Hz (Hertz), s (seconds), m (meters), dB (deziBel)

clear, clc
% communication model parameters
N0      = -170;          % in dBm/Hz
B       = 1e06;          % Hz. channel bandwidth
alpha0  = -50;           % channel power gain in dB at the reference distanc 1m
sig0    = sqrt(B*N0);    % unitless. Noise power at the receiver 

% radar sensing model 
Gp      = 0.1*B;        % signal processing gain
beta0   = -47;      % channel power gain in dB at the reference distanc 1m

% other parameters
H      = 200;    % in m. Flight height of the qudcopter, which we assum to be constant.
Vstr   = 20;  % in m/s. Fixed flight speed for the initialization process.
Lx     = 1500;  % in m. Length of the x dimension
Ly     = 1500;  % in m. Length of the y dimension 
Tf     = 1.5;   % in s. flight time? (used for the energy function).
Th     = 1;  % in s. Hover time? (used for the energy function).
Etotal = 35e03;
Vmax   = 30;  % in m/s. Maximum flight speed.
eta    = .8; % unitless. Weighting factor to weigh between the sensing and communication objective. Higher values favor sensing, lower values communication.  

% trajecotry parameters
M    = 1;      % amount of stages 
Nstg = 25;  % amount of points in each stage (the same for all stages except for the last one 
mu   = 5;     % hover steps (hover every mu steps)
Kstg = floor(Nstg/mu); % amount of hover points in each stage (the same for all stages except for the last one)
