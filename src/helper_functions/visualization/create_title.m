function ttl = create_title(varargin)
%CREATE_TITLE Summary of this function goes here
%   Detailed explanation goes here

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
    
    str_ttl = "(eta-" + eta + ")(w_{star}-" + w_star + ")(N_{stg}-" + N_stg + ")(mu-" + mu+")";  
    ttl = title(axs,str_ttl);
end

end
