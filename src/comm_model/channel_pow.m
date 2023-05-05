function h = channel_pow(alpha_0, dc)
%CHAMNEL_POW function for the channelpower of the communication channel
%   alpha_0 == channel power at referenve distance 1 m
%   dc == distance to the communication user

h = alpha_0./(dc^2);
end

