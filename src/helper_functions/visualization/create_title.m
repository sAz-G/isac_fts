%------------------------------------------------------------------------
% FUNCTION NAME: create_title
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Creates a title for the trajectory map.
%
% INPUTS:
% varargin - Contains parameters and axis.
%
% OUTPUTS:
% ttl - The title.
%
% USAGE: ttl = create_title(varargin)
%-----------------------------------------------------------------------

function ttl = create_title(varargin)
    def_ttl = "Default";

    if nargin == 0 || nargin == 1
        axs = varargin{1};
        ttl = title(axs, def_ttl);
    elseif nargin == 2
        axs = varargin{1}; % Axis
        prms = varargin{2}; % Parameters

        eta = prms.sim.eta;
        w_star = prms.sim.w_star;
        N_stg = prms.sim.N_stg;
        mu = prms.sim.mu;
        iter = prms.sim.iter;
        str_ttl = "Trajectory: ($\eta=" + eta + ")(w=" + w_star + ")(N_{m}=" + N_stg + ")(\mu=" + mu + ")(" + "iter=" + iter + ")$";  
        ttl = title(axs, str_ttl, 'Interpreter', 'latex');
    end
end
