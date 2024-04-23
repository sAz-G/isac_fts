%------------------------------------------------------------------------
% FUNCTION NAME: plot_map
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Creates a trajectory map.
%
% INPUTS:
% S       - Trajectory.
% s_b     - Charging station.
% s_t     - Sensing target.
% s_t_est - Estimated positions of sensing target.
% s_c     - Communication user.
% varargin - Further parameters.
%
% OUTPUTS:
% p_s     - Plot handle for trajectory.
% p_h     - Plot handle for hovering.
% p_b     - Plot handle for charging station.
% p_t     - Plot handle for sensing target.
% p_c     - Plot handle for communication user.
% p_t_hat - Plot handle for estimated target position.
% fig     - Figure handle.
%
% USAGE: [p_s, p_h, p_b, p_t, p_c, p_t_hat, fig] = plot_map(S, s_b, s_t, s_t_est, s_c, varargin)
%-----------------------------------------------------------------------

function [p_s, p_h, p_b, p_t, p_c, p_t_hat, fig] = plot_map(S, s_b, s_t, s_t_est, s_c, varargin)
    params = [];
    if ~isempty(varargin)
        params = varargin{1};
    end
    
    fig = figure;
    xlabel("x [m]")
    ylabel("y [m]")
    xlim([0 1500])
    ylim([0 1500])
    hold on
    grid on
    
    % Plot trajectory
    p_s = plot(S(1,:), S(2,:), 'bo-', 'MarkerSize', 2);
    
    % Plot charging station
    p_b = plot(s_b(1), s_b(2), 'g*', 'MarkerSize', 10);
    
    % Plot sensing target
    p_t = plot(s_t(1), s_t(2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % Plot estimated target position
    p_t_hat = plot(s_t_est(1,:), s_t_est(2,:), 'b-', 'LineWidth', 1.1, 'Marker', 'diamond', 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none');
    
    % Plot communication user
    p_c = plot(s_c(1), s_c(2), 'cs', 'MarkerSize', 10, 'MarkerFaceColor', 'c');
    
    % Set legend
    legend({'S', 's_b', 's_t', 's_{est}', 's_c'}, 'Location', 'southeast');
    
    % Set axis ticks
    set(gca, 'xtick', 0:250:1500)
    set(gca, 'ytick', 0:250:1500)

    if ~isempty(params)
        mu = params.sim.mu;
        
        % Plot hovering
        p_h = plot(S(1,mu:mu:end), S(2,mu:mu:end), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
        legend({'S', 's_b', 's_t', 's_{est}', 's_c', 'S_h'}, 'Location', 'southeast');
        
        % Set title
        create_title(gca, params);
    else
        % Set title
        create_title(gca);
    end
end
