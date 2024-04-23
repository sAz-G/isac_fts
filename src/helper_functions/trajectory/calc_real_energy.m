%------------------------------------------------------------------------
% FUNCTION NAME: calc_constraint_energy
% AUTHOR: Sharif Azem     (TU-Darmstadt department 18, sAz-G on github)
%         Markus Krantzik (TU-Darmstadt department 18, mardank on github)
%
% DESCRIPTION:  calculates the used energy, which is needed for the
% multi stage problem problem
%
% INPUTS:
%  S - trajectory 
%  s_s - last point of the previous stage
%  params - predefined parameters
%
% OUTPUTS:
%       E_used - constraint energy 
%
% USAGE: E_used = calc_real_energy(S, s_s, params)
%-----------------------------------------------------------------------

function E_used = calc_real_energy(S, s_s, params)    
    % enrgy parameters 
    s     = params.energy.s; 
    A     = params.energy.A;
    rho   = params.energy.rho;
    D_0   = params.energy.D_0;
    v_0   = params.energy.v_0;
    U_tip = params.energy.U_tip;
    P_0   = params.energy.P_0;
    P_I   = params.energy.P_I;
    
    T_f = params.sim.T_f; 
    T_h = params.sim.T_h; 
    K_stg = params.sim.K_stg;
    
    V = calc_velocity(S, s_s, params);
    V_norm = norms(V, 2, 1);
    
    % power and energy 
    flight_power = power_model(V_norm);
    hover_power  = power_model(zeros(1, K_stg));
    E_used       = T_f * sum(flight_power) + T_h * sum(hover_power);
    
    % a function to calculate the power 
    function P = power_model(V_norm)
        P = P_0 * (1 + 3*V_norm.^2/U_tip^2) + P_I * sqrt(sqrt(1 + V_norm.^4/(4*v_0^4)) - V_norm.^2/(2*v_0^2)) + 1/2 * D_0 * rho * s * A * V_norm.^2;
    end
end

