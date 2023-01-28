%% predefined hyper parameters for the algorithm
% below the used units are provided. Using known already units.
% Hz (Hertz), s (seconds), m (meters), dB (deziBel)

%% Communication model parameters
N_0     = db2pow(-170) * 1e-3;  % [W/Hz]
B       = 1e06;                 % channel bandwidth [Hz]
alpha_0 = db2pow(-50);          % channel power gain in dB at the reference distanc d_c(n) = 1m
sigma_0 = sqrt(B * N_0);        % unitless. Noise power at the receiver 
P       = db2pow(20) * 1e-3;    % transmit power [dB]

% radar sensing model 
G_p      = 0.1 * B;     % signal processing gain [Hz]
beta_0   = db2pow(-47); % channel power gain in dB at the reference distanc d_s(k) = 1m [W]

%% Energy parameters
% Parameters for the derivation of O_I and P_0
% delta_eng = 0.012;
% W     = 20;       % aircraft weight in Newton
% A     = 0.503;    % rotor_disc area in [m^2]
% Sigma = 300;  % blade angular velocity in [radians/second]
% R     = 0.4;      % rotor radius in [m]
% k     = 0.1;      % incremental correction factor to include power;

E_total = 35e03;    % total available energy [J]
P_0     = 80;       % blade profile power [W] (P_0 = delta_eng/8 * rho * s * A * Sigma^3 * R^3)
P_I     = 88.6;     % induced power in hovering status [W] (P_I = (1 + k) * (W^(3/2))/sqrt(2*rho * A))
U_tip   = 120;      % tip speed of the rotor [m/s]
v_0     = 4.03;     % mean rotor induced velocity in forward flying [m/s]
D_0     = 0.6;      % Fuselage drag ratio
s       = 0.05;     % rotor solidity [m^3]
rho     = 1.225;    % air density [kg/m^3]
A       = 0.503;    % rotor_disc area in [m^2]

%% other parameters
H       = 200;    % flight height of the qudcopter, which we assum to be constan [m]
V_str   = 20;    % fixed flight speed for the initialization process [m/s]
L_x     = 1500;  % length of the x dimension [m]
L_y     = 1500;  % length of the y dimension [m]
T_f     = 1.5;   % flying duration/ flight time [s] (used for the energy function)
T_h     = 1;  % hovering duration/ hover time [s] (used for the energy function)
V_max   = 30;  % maximum flight speed [m/s]
a = 10; % pre-determined constant related to the system setting

% Optimization parameter
% eta     = .8; % unitless. Weighting factor to weigh between the sensing and communication objective. Higher values favor sensing, lower values communication.

%% trajecotry parameters
% M    = 1;      % amount of stages 
% Nstg = 25;  % amount of points in each stage (the same for all stages except for the last one 
% mu   = 5;     % hover steps (hover every mu steps)
% Kstg = floor(Nstg/mu); % amount of hover points in each stage (the same for all stages except for the last one)
