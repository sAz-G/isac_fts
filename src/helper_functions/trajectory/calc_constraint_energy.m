%------------------------------------------------------------------------
% FUNCTION NAME: calc_constraint_energy
% AUTHOR: Sharif Azem (sAz-G on github), Markus Krantzik (mardank on github)
%
%
% DESCRIPTION:  calculates the constraint energy, which is needed for the
% optimization problem
%
% INPUTS:
%  V - velocity 
%  delta - opt. variable 
%  params - predefined parameters
%
% OUTPUTS:
%       E_const - constraint energy 
%
% USAGE: E_const = calc_constraint_energy(V,delta,params)
%-----------------------------------------------------------------------


function E_const = calc_constraint_energy(V,delta,params)

% energy parameters 
s     = params.energy.s; 
A     = params.energy.A;
rho   = params.energy.rho;
D_0   = params.energy.D_0;
U_tip = params.energy.U_tip;
P_0   = params.energy.P_0;
P_I   = params.energy.P_I;

T_f   = params.sim.T_f;
T_h   = params.sim.T_h;
K_stg = params.sim.K_stg;

% power 
power_sum1 = P_0.*(size(V,2) + 3./(U_tip.^2).*sum( pow_pos( norms(V,2,1),2) ) );
power_sum1 = power_sum1 + 0.5*D_0*rho*s*A*sum(pow_pos(norms(V,2,1),3));
power_sum2 = P_I.*sum(delta);
power_sum3 = K_stg.*(P_0+P_I);

% energy
E_const     = T_f*(power_sum1+power_sum2)+T_h*power_sum3;
 

end