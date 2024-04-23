%------------------------------------------------------------------------
% FUNCTION NAME: create_title
% AUTHOR: Sharif Azem     (TU-Darmstadt department 18, sAz-G on github)
%         Markus Krantzik (TU-Darmstadt department 18, mardank on github)
%
% DESCRIPTION:  creates a title for the trajectory map
%
% INPUTS:
% varargin - contains parameters and axis
%
% OUTPUTS:
%       ttl - the title 
%
% USAGE:  ttl = create_title(varargin)
%-----------------------------------------------------------------------

function ttl = create_title(varargin)

def_ttl = "Default";
if nargin == 0 % enter here just in case the user does not insert a title 
    ttl = title(def_ttl);
elseif nargin == 1 % enter here just in case the user does not insert a title 
    axs = varargin{1};
    ttl = title(axs , def_ttl);
elseif nargin == 2 % create a title with the given parameters
    axs = varargin{1}; % axis 
    prms = varargin{2}; % parameters

    eta = prms.sim.eta;
    w_star = prms.sim.w_star;
    N_stg = prms.sim.N_stg;
    mu = prms.sim.mu;
    iter = prms.sim.iter;
    str_ttl = "Trajectory: ($\eta=" + eta + ")(w=" + w_star + ")(N_{m}=" + N_stg + ")(\mu=" + mu+")(" + "iter=" + iter + ")$";  
    ttl = title(axs,str_ttl, 'Interpreter', 'latex');
    
end

end

