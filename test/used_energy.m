function out_energy = used_energy(V, K_m)
%USED_ENERGY Summary of this function goes here
%   Detailed explanation goes here
    
    Em_1 = sum(P_0 * (1 + 3* norm(V).^2/(U_tip.^2)));
    Em_1 = Em_1 + P_I * sqrt(sqrt(1 + norm(v).^4/(4*v_0^4)) - norm(V)^2./(2*v_0^2));
    Em_1 = Em_1 + 0.5 * D_0 * rho * s * A * norm(V)^3;

    out_energy = T_f * sum(Em_1) + T_h * K_m * (P_0 + P_I);
end

