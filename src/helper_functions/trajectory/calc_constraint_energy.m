%------------------------------------------------------------------------
% FUNCTION NAME: calc_constraint_energy
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Calculates the constraint energy for the optimization problem.
%
% INPUTS:
%   V       - Velocity matrix.
%   delta   - Optimal variable.
%   params  - Predefined parameters.
%
% OUTPUTS:
%   E_const - Constraint energy.
%
% USAGE: E_const = calc_constraint_energy(V, delta, params)
%-----------------------------------------------------------------------

function E_const = calc_constraint_energy(V, delta, params)
    % Extract energy parameters from params
    s     = params.energy.s; 
    A     = params.energy.A;
    rho   = params.energy.rho;
    D_0   = params.energy.D_0;
    U_tip = params.energy.U_tip;
    P_0   = params.energy.P_0;
    P_I   = params.energy.P_I;

    % Extract simulation parameters from params
    T_f   = params.sim.T_f;
    T_h   = params.sim.T_h;
    K_stg = params.sim.K_stg;

    % Calculate power sums
    power_sum1 = P_0 * (size(V, 2) + 3 / (U_tip^2) * sum(norms(V, 2, 1).^2));
    power_sum1 = power_sum1 + 0.5 * D_0 * rho * s * A * sum(norms(V, 2, 1).^3);
    power_sum2 = P_I * sum(delta);
    power_sum3 = K_stg * (P_0 + P_I);

    % Calculate constraint energy
    E_const = T_f * (power_sum1 + power_sum2) + T_h * power_sum3;
end
