function E_used = calc_real_energy(K_stg, S, S_s, params)
    % calc_energy calculates the energy, given the previous energy, and the
    % energy used in each stage, which is predefined.
    % Arguments: 
    % E_min is the energy used in each stage excluding the final
    % stage.
    % E_prev is the previous energy state.
    % Returns: 
    % E_new is the new energy state
    
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

    V_froms_S = calc_velocity_from_trajectory(S, S_s, params);
    V_abs = norms(V_froms_S, 2, 1);
    
    Energy_abs_vec = P_enery(V_abs);
    Energy_zero_vec = P_enery(zeros(1, K_stg));
    E_used = T_f * sum(Energy_abs_vec) + T_h * sum(Energy_zero_vec);

    function out_P = P_enery(V_abs_in)
        out_P = P_0 * (1 + 3*V_abs_in.^2/U_tip^2) + P_I * sqrt(sqrt(1 + V_abs_in.^4/(4*v_0^4)) - V_abs_in.^2/(2*v_0^2)) + 1/2 * D_0 * rho * s * A * V_abs_in.^2;
    end
end

