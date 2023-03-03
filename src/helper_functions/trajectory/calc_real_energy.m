function E_used = calc_real_energy(S, s_s, params)
    % calc_energy calculates the energy, given the previous energy, and the
    % energy used in each stage, which is predefined.
    % Arguments: 
    % E_min is the energy used in each stage excluding the final
    % stage.
    % E_prev is the previous energy state.
    % Returns: 
    % E_new is the new energy state
    
    % enrgy parameters 
    s     = params.energy.s; 
    A     = params.energy.A;
    rho   = params.energy.rho;
    D_0   = params.energy.D_0;
    v_0   = params.energy.v_0;
    U_tip = params.energy.U_tip;
    P_0   = params.energy.P_0;
    P_I   = params.energy.P_I;
    
    T_f = params.sim.T_f; % flight time
    T_h = params.sim.T_h; % hover time 
    K_stg = params.sim.K_stg;
    
    V = calc_velocity(S, s_s, params);
    V_norm = norms(V, 2, 1);
    
    flight_power = power_model(V_norm);
    hover_power  = power_model(zeros(1, K_stg));
    E_used       = T_f * sum(flight_power) + T_h * sum(hover_power);

    function P = power_model(V_norm)
        P = P_0 * (1 + 3*V_norm.^2/U_tip^2) + P_I * sqrt(sqrt(1 + V_norm.^4/(4*v_0^4)) - V_norm.^2/(2*v_0^2)) + 1/2 * D_0 * rho * s * A * V_norm.^2;
    end
end

