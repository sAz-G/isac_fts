function E_used = calc_constraint_energy(K_stg,V,delta,params)
% calc_constraint_energy calculates Calculate the energy for a given stage 
% This function uses the convex energy function, which approximates the
% energy in the constraints
%   K_stg == number of hovering points
%   S == trajectory of the flight
%   S_s == start point of the trajectory
%   params == struct with simulation parameters

s     = params.energy.s; 
A     = params.energy.A;
rho   = params.energy.rho;
D_0   = params.energy.D_0;
U_tip = params.energy.U_tip;
P_0   = params.energy.P_0;
P_I   = params.energy.P_I;

T_f = params.sim.T_f;
T_h = params.sim.T_h;

Em_sum1 = sum(P_0 * (1 + 3 * pow_pos(norms(V,2,1), 2)/(U_tip.^2)) + 0.5 * D_0 * rho * s * A * pow_pos(norms(V, 2, 1), 3));
Em_sum2 = sum(P_I*delta);
Em_sum3 = K_stg* (P_0 + P_I);
E_used = T_f * Em_sum1 + T_f * Em_sum2 + T_h * Em_sum3;
 

end

