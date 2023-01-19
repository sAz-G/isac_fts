function [p_s, p_b, p_t, p_c, p_t_hat] = plot_map(S, s_b, s_t, s_t_est, s_c)
%PLOT_MAP plots the map for a given trajectory and base station, target and
%communication user position.
% arguments: 
% S    - quad trajectory.
% s_b  - base station position.
% s_t  - target position. 
% s_c - communication user position.

% map 
    figure
    hold on
    grid on
    
    % plot trajectory 
    p_s = plot(S(1,:), S(2,:));
    p_s.LineStyle = 'none';
    p_s.Marker = 'o';
    p_s.LineStyle = 'none';
    p_s.MarkerEdgeColor = 'b';
    
    % plot base station 
    p_b = plot(s_b(1), s_b(2));
    p_b.LineStyle = 'none';
    p_b.Marker = '*';
    p_b.MarkerEdgeColor = 'g';
    
    % plot target position 
    p_t = plot(s_t(1), s_t(2));
    p_t.LineStyle = 'none';
    p_t.Marker = '^';
    p_t.MarkerSize = 10;
    p_t.MarkerFaceColor = 'r';
    p_t.MarkerEdgeColor = 'r';
    
    % plot estimated target position 
    p_t_hat = plot(s_t_est(1), s_t_est(2));
    p_t_hat.LineStyle = 'none';
    p_t_hat.LineWidth = 1.1;
    p_t_hat.Marker = 'diamond';
    p_t_hat.MarkerSize = 10;
    p_t_hat.MarkerFaceColor = 'none';
    p_t_hat.MarkerEdgeColor = 'b';
    
    % plot communication user 
    p_c = plot(s_c(1), s_c(2));
    p_c.LineStyle = 'none';
    p_c.Marker = 'square';
    p_c.MarkerSize = 10;
    p_c.MarkerFaceColor = 'c';
    p_c.MarkerEdgeColor = 'c';
    
    legend('S', 's_b', 's_t','s_{est}', 's_c')
end

