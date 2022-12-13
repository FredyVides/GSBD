function [Ad, Bd, Q, R] = HVACcreator(nz, tol, Ts)

C = 1/1.4;
Ca=C*1.1;
fac=1;
Cf = C*fac;

% nz = 4;

[Adj, nw] = HVAC_Adj(nz, C, Cf, Ca);
L = diag(sum(Adj,2))-Adj;

% C_z = 250000 ; %Capacitance J/Kg
% C_w = 4660000;

C_z = 250 ; %Capacitance J/Kg
C_w = 4660;

na=nz+nw;

R_z = eye(nz).*1/C_z ;  %Resistance of walls
R_w = eye(nw).*1/C_w;

%Use algorithm 1 from paper

Lgx=L(1:na,1:na);
Lga=L(1:nw,end);

A_cont = -diag([diag(R_w);diag(R_z)])*Lgx;
Aww = A_cont( 1:nw,  1:nw);
Azz = A_cont(nw+1:nw+nz, nw+1:nw+nz);
Azw = A_cont(nw+1:nw+nz,  1:nw);
Awz = A_cont( 1:nw, nw+1:nw+nz);
A = [Azz, Azw; Awz, Aww];
%A_cont = ones(13)
B_a_cont = -R_w*Lga;
B_z_cont = R_z;

% This is if we dont ignore consttnt ambient temp and 0 disruption event
% B_cont = [zeros([size(B_a_cont,1),size(B_z_cont,2)]),B_a_cont,zeros([size(B_a_cont,1),size(B_z_cont,2)]);B_z_cont, zeros([size(B_z_cont,1),size(B_a_cont,2)]), B_z_cont];


% This is if we ignore constant ambient temp and 0 disruption event
B_cont = [zeros([nw,size(B_z_cont,2)]);B_z_cont];

m = size(B_cont,2);
n = size(A_cont,2);
%
B = [B_z_cont; zeros([nw,size(B_z_cont,2)])];
B=[B,zeros(nw+nz,nw+nz-size(B,2))];

C_cont = diag([zeros([1,nw]),ones([1,nz])]);
C = diag([ones([1,nz]), zeros([1,nw])]);

D_cont = zeros(size(B_cont));
D = zeros(nw+nz);

%Cont system
% sys = ss(A_cont,B_cont,C_cont,D_cont);

sys = ss(A,B,C,D);

%Discrete System
sysd = c2d(sys,Ts);
Ad = sysd.A; Ad(abs(Ad)<tol) = 0;
Bd = sysd.B; Bd(abs(Bd)<tol) = 0;



Q = diag([ones([nz,1]);0.01*ones([nw,1])]);
R = blkdiag(0.1*eye(m), eye(n-m));
end