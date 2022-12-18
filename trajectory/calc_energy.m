function [E_new] = calc_energy(E_prev, E_min)
% calc_energy calculates the energy, given the previous energy, and the
% energy used in each stage, which is predefined.
% Arguments: 
% E_min is the energy used in each stage excluding the final
% stage.
% E_prev is the previous energy state.
% Returns: 
% E_new is the new energy state

E_new = E_prev - E_min;

end

