function bound = rate_boundary(s_c,s_start, N_traj, params,tp)
%RATE_BOUNDARY Summary calculates the minimum or maximum achievable rate
%values This can be used to bound the values of the objective function,
%which includes an approximation of the rate using its taylor series
%expansion.
%   
%Input: 
% s_c       - position of the communicatin user
% s_start   - current position of the quadcopter
% N_traj    - the total amount of points of the whole trajectory.
% params    - simulation parameters and constants
% tp        - type of boundary, minimum or maximum
% output:
% bound     - the maximum or minimum value

T_f = params.sim.T_f; % flight time between points 
N_stg = params.sim.N_stg; % future trajectory points (N_stg <= N_traj), (S_stg is a subset of S_traj)

if strcmp(tp, 'max') % calculate maximum value
    S_bound = zeros(2,params.sim.N_stg);
    % calculat the trajectory that gives the maximum avveraged data rate 
    max_stp = params.sim.V_max.*params.sim.T_f; % maximum distance between two trajectory points
    user_dist = norm(s_c-s_start,2); % distance to the communication user
    N_dist = floor(user_dist./max_stp); % amount of points to the communication user, when flying with the maximum speed
    
    if N_dist >= N_stg % if the quad reaches the communication user before flying N_stg points
        S_bound(:,1:N_stg) = s_start + params.sim.V_max.*(1:N_stg).*(s_c - s_start)./norm(s_c - s_start)*T_f;
    elseif N_dist < N_stg % if the quad can not reach the communication user flying N_stg points 
        % calculate a trajectory to the user
        S_bound(:,1:N_dist) = s_start + params.sim.V_max.*(1:N_dist).*(s_c - s_start)./norm(s_c - s_start)*T_f;
        % after reaching the user, all points are at the same position as
        % the communication user
        S_bound(:,N_dist+1:N_stg) = s_c.*ones(2,abs(N_dist-N_stg));
    end
elseif strcmp(tp, 'min') % calculate minimum value
    % the trajectory for the minimum value always points away from the
    % communication user. The quad flies with V_max to get as far as
    % possible from the communication user
    S_bound(:,1:N_stg) = s_start-params.sim.V_max .*(1:N_stg).* (s_c - s_start)./norm(s_c - s_start)*T_f;
end

bound = avg_data_rate(S_bound,s_c,params, N_traj); % claculate the boundary value
end

