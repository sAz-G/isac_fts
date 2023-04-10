function [S_init, V_init] = init_trajectory(s_start, s_end,N, params)
%INIT_TRAJECTORY computes the initial trajectory for the iterative process
% arguments:
% s_start - starting point, the point at which the quad starts.
% N       - the amount of trajectory points at the mth stage.
% s_end   - the s_end point of the trajectory between the communication user and the
%           estimated sensing target (the estimated position obtained from the
% radar).
% V_str   - the speed of the quad. It is constant throughout the
% initial trajectory.
% T_f     - flight time between points.
% output:
% S_init  - initial trajectory.
% V_out   - velocity throughout the trajectory.

V_str = params.sim.V_str;
T_f   = params.sim.T_f;
L_x = params.sim.L_x;
L_y = params.sim.L_y;

S_init = zeros(2,N);
V_init  = ones(2,N).*V_str.*(s_end - s_start)./norm(s_end - s_start);

S_init(:,1:N) = s_start + V_str.*(1:N).*(s_end - s_start)./norm(s_end - s_start)*T_f;

% check if the calculated trajectory is valid. if not, calculate a new one.
s_e = S_init(:,end); % get last position at the initial trajectory


if (s_e(1) > L_x)&&(s_e(2) > L_y)                                    % x position x out of bounds
    
    if (abs(s_e-[L_x; s_e(2)]) > abs(s_e-[s_e(1); L_y]))
        x_new  = L_x;                                                    % calc new x
        s_e(1) = x_new;                             % calc new x
        grad   = s_start - s_end;                   % calc diff to obtain slope
        grad   = grad(2)./grad(1);                  % calc slop
        y_new  = s_e(2) + grad.*(x_new-s_e(1));   % calc y (no need to check boundary since x boundary is checked, and modification is linear)
        s_e(2) = y_new;                             % set new x
    else
        y_new     = L_y;
        s_e(2)    = y_new;
        grad_inv  = s_start - s_end;
        grad_inv  = grad_inv(1)./grad_inv(2);
        x_new     = s_end(1) + grad_inv.*(y_new-s_e(2));
        s_e(1)    = x_new;
    end
    
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N); % obtain linear trajectory. Note velocity is different than input
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif (s_e(1) > L_x)&&(s_e(2) < 0)
    if (abs(s_e-[L_x; s_e(2)]) > abs(s_e-[s_e(1); 0]))
        x_new  = L_x;                                                    % calc new x
        s_e(1) = x_new;                             % calc new x
        grad   = s_start - s_end;                   % calc diff to obtain slope
        grad   = grad(2)./grad(1);                  % calc slop
        y_new  = s_e(2) + grad.*(x_new-s_e(1));   % calc y (no need to check boundary since x boundary is checked, and modification is linear)
        s_e(2) = y_new;                             % set new x
    else
        y_new     = 0;
        s_e(2)    = y_new;
        grad_inv  = s_start - s_end;
        grad_inv  = grad_inv(1)./grad_inv(2);
        x_new     = s_end(1) + grad_inv.*(y_new-s_e(2));
        s_e(1)    = x_new;
    end
    
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N); % obtain linear trajectory. Note velocity is different than input
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif (s_e(1) < 0)&&(s_e(2) > L_y)
    if (abs(s_e-[0; s_e(2)]) > abs(s_e-[s_e(1); L_y]))
        x_new  = 0;                                                    % calc new x
        s_e(1) = x_new;                             % calc new x
        grad   = s_start - s_end;                   % calc diff to obtain slope
        grad   = grad(2)./grad(1);                  % calc slop
        y_new  = s_e(2) + grad.*(x_new-s_e(1));   % calc y (no need to check boundary since x boundary is checked, and modification is linear)
        s_e(2) = y_new;                             % set new x
    else
        y_new     = L_y;
        s_e(2)    = y_new;
        grad_inv  = s_start - s_end;
        grad_inv  = grad_inv(1)./grad_inv(2);
        x_new     = s_end(1) + grad_inv.*(y_new-s_e(2));
        s_e(1)    = x_new;
    end
    
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N); % obtain linear trajectory. Note velocity is different than input
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif (s_e(1)  < 0)&&(s_e(2)  < 0)
    if (abs(s_e-[0; s_e(2)]) > abs(s_e-[s_e(1); 0]))
        x_new  = 0;                                                    % calc new x
        s_e(1) = x_new;                             % calc new x
        grad   = s_start - s_end;                   % calc diff to obtain slope
        grad   = grad(2)./grad(1);                  % calc slop
        y_new  = s_e(2) + grad.*(x_new-s_e(1));   % calc y (no need to check boundary since x boundary is checked, and modification is linear)
        s_e(2) = y_new;                             % set new x
    else
        y_new     = 0;
        s_e(2)    = y_new;
        grad_inv  = s_start - s_end;
        grad_inv  = grad_inv(1)./grad_inv(2);
        x_new     = s_end(1) + grad_inv.*(y_new-s_e(2));
        s_e(1)    = x_new;
    end
    
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N); % obtain linear trajectory. Note velocity is different than input
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif s_e(1) > L_x
    x_new  = L_x;                               % calc new x
    s_e(1) = x_new;                             % calc new x
    grad   = s_start - s_end;                   % calc diff to obtain slope
    grad   = grad(2)./grad(1);                  % calc slop
    y_new  = s_end(2) + grad.*(x_new-s_e(1));   % calc y (no need to check boundary since x boundary is checked, and modification is linear)
    s_e(2) = y_new;                             % set new x
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N); % obtain linear trajectory. Note velocity is different than input
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif s_e(1) < 0
    x_new = 0;
    s_e(1) = x_new;
    grad   = s_start - s_end;
    grad   = grad(2)./grad(1);
    y_new  = s_end(2) + grad.*(x_new-s_e(1));
    s_e(2) = y_new;
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N);
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif s_e(2) > L_y
    y_new     = L_y;
    s_e(2)    = y_new;
    grad_inv  = s_start - s_end;
    grad_inv  = grad_inv(1)./grad_inv(2);
    x_new     = s_end(1) + grad_inv.*(y_new-s_e(2));
    s_e(1)    = x_new;
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N);
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif s_e(2) < 0
    y_new    = 0;
    s_e(2)   = y_new;
    grad_inv = s_start - s_end;
    grad_inv = grad_inv(1)./grad_inv(2);
    x_new    = s_end(1) + grad_inv.*(y_new-s_e(2));
    s_e(1)   = x_new;
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N);
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
end

end

