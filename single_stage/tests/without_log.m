clear
close all

%% Simulation parameter

% log(1+1/x) = log((1+x)/x) = -log(x/(1+x)) = -log(1 - 1/(1+x)) = -log(1 - inv_pos(1+x))

ub = 1e7

cvx_begin
    cvx_solver mosek
    variable x;
    maximize(1 - inv_pos(1+x)); % optimize without log
    subject to
        x >= 0; % exponatial cone
        ub >= x;
cvx_end

x_with_sub = x;

cvx_begin
    cvx_solver mosek
    variable x;
    maximize(log(1 - inv_pos(1+x))); % optimize without log
    subject to
        x >= 0; % exponatial cone
        ub >= x;
cvx_end

x_without_sub = x;

x_with_sub
x_without_sub



