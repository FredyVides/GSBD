clc
clear
close all
%% HVAC loading
% either 4, 9, 16, 25, 36
zone = 4; 
data = load(['zone', num2str(zone),'.mat']);
Ap = data.Ad; Ap(abs(Ap)<1e-6) = 0;
Bp = data.Bd; Bp(abs(Ap)<1e-6) = 0;
Qp = data.Q;
Rp = data.R; 

%% generating constriant matrices
% length of original vector u (control input vector)
nc = size(data.B_cont,2); 

% length of x (state vector)
n = length(Ap); 

% set 2*nc (all original control inputs are bounded, i.e., a <= u_i <=b)
% or set 2*n (all states are bounded, i.e., a <= x_i <=b)
num_constraints = 2*nc;
Gp = zeros(num_constraints, length(Ap)); 
for i = 2:2:num_constraints
    Gp(i-1,i/2) = 1; % meaning: 1*x_i <= b
    Gp(i,i/2) = -1; % meaning: -1*x_i <= -a or equivalently a <= x_i
end

% Then we use sbd:  T = sbd(Ap, Bp, Qp, Rp, ....) without further
% perturbing the matrices
