%------------------------------------------------------------------------
% FUNCTION NAME: plot_map
% AUTHOR: Sharif Azem     (TU-Darmstadt department 18, sAz-G on github)
%         Markus Krantzik (TU-Darmstadt department 18, mardank on github)
%
% DESCRIPTION: creates a trajectory map
%
% INPUTS:
%  S   - trajectory 
%  s_b - charging station 
%  s_t - sensing target
%  s_t_est - estimated positions of sensing target
%  s_c - communicatoin user
%  varargin - further parameters
%
% OUTPUTS:
%       the plots and the figure 
%
% USAGE: [p_s,p_h, p_b, p_t, p_c, p_t_hat,fig] = plot_map(S, s_b, s_t, s_t_est, s_c, varargin)
%-----------------------------------------------------------------------

function [p_s,p_h, p_b, p_t, p_c, p_t_hat,fig] = plot_map(S, s_b, s_t, s_t_est, s_c, varargin)
    
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
    
    % plot trajectory 
    p_s = plot(S(1,:), S(2,:));
    p_s.Marker = 'o';
    p_s.MarkerEdgeColor = 'b';
    p_s.MarkerEdgeColor = 'b';
    p_s.MarkerSize      = 2;

    % plot charging station 
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
    p_t_hat = plot(s_t_est(1,:), s_t_est(2,:));
    p_t_hat.LineStyle = '-';
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
    
    l = legend();
    l.Location = 'southeast';
    
    set(gca,'xtick',[0:250:1500])
    set(gca,'ytick',[0:250:1500])

    if isempty(params)
         l.String = {'S','s_b', 's_t','s_{est}', 's_c'};         
        create_title(gca);
    else
        % plot hovering
        mu = params.sim.mu;
        p_h = plot(S(1,mu:mu:end), S(2,mu:mu:end));
        p_h.LineStyle = 'none';
        p_h.Marker = 'o';
        p_h.LineStyle = 'none';
        p_h.MarkerEdgeColor = 'b';    
        p_h.MarkerFaceColor = 'k';
        p_h.MarkerSize      = 4;

        l.String = {'S','s_b','s_t', 's_{est}','s_c', 'S_h'};
        create_title(gca, params);
    end
    
    
end

