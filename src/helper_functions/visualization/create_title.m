%------------------------------------------------------------------------
% FUNCTION NAME: create_title
% AUTHOR: Sharif Azem
%         Markus Krantzik
%
% DESCRIPTION: creates the title of the trajectory map, which is created
% using the function plot_map
%
% INPUTS:
%   varargin - current axis and a struct 
%
% OUTPUTS:
%   ttl - the title
%
% USAGE:  ttl = create_title(gca,params);
%
%------------------------------------------------------------------------

function ttl = create_title(varargin)

def_ttl = "Default";
if nargin == 0
    ttl = title(def_ttl);
elseif nargin == 1
    axs = varargin{1};
    ttl = title(axs , def_ttl);
elseif nargin == 2
    axs = varargin{1};
    prms = varargin{2};

    eta = prms.sim.eta;
    w_star = prms.sim.w_star;
    N_stg = prms.sim.N_stg;
    mu = prms.sim.mu;
    iter = prms.sim.iter;
    str_ttl = "Trajectory: ($\eta=" + eta + ")(w=" + w_star + ")(N_{m}=" + N_stg + ")(\mu=" + mu+")(" + "iter=" + iter + ")$";  
    ttl = title(axs,str_ttl, 'Interpreter', 'latex');
    
end

end

