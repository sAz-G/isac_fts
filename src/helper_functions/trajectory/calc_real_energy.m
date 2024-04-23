%------------------------------------------------------------------------
% FUNCTION NAME: calc_real_energy
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Calculates the used energy for the multi-stage problem.
%
% INPUTS:
%   S       - Trajectory matrix.
%   s_s     - Last point of the previous stage.
%   params  - Predefined parameters.
%
% OUTPUTS:
%   E_used  - Total energy used.
%
% USAGE: E_used = calc_real_energy(S, s_s, params)
%-----------------------------------------------------------------------

function E_used = calc_real_energy(S, s_s, params)    
    % Extract energy parameters from params
    s     = params.energy.s; 
    A     = params.energy.A;
    rho   = params.energy.rho;
    D_0   = params.energy.D_0;
    v_0   = params.energy.v_0;
    U_tip = params.energy.U_tip;
    P_0   = params.energy.P_0;
    P_I   = params.energy.P_I;
    
    % Extract simulation parameters from params
    T_f   = params.sim.T_f; 
    T_h   = params.sim.T_h; 
    K_stg = params.sim.K_stg;
    
    % Calculate velocity
    V = calc_velocity(S, s_s, params);
    V_norm = norms(V, 2, 1);
    
    % Calculate flight and hover power
    flight_power = power_model(V_norm);
    hover_power  = power_model(zeros(1, K_stg));
    
    % Calculate total energy used
    E_used = T_f * sum(flight_power) + T_h * sum(hover_power);
    
    % Nested function to calculate power
    function P = power_model(V_norm)
        P = P_0 * (1 + 3 * V_norm.^2 / U_tip^2) + P_I * sqrt(sqrt(1 + V_norm.^4 / (4 * v_0^4)) - V_norm.^2 / (2 * v_0^2)) + 1/2 * D_0 * rho * s * A * V_norm.^2;
    end
end
