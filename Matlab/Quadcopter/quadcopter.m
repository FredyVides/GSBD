clc
clear
close all

%% constants
g = 9.81; % gravitational constant
mass = 0.65; % mass of quadcopter
Ixx = 7.5e-3; % moment of inersia
Iyy = Ixx; % moment of inersia
Izz = 1.3e-2; % moment of inersia

n = 12; % number of states
m = 4; % number of inputs
%% matrices from quadcopter

Nc = 8; % number of quadcopters (variable)

Ahat0 = [
    0 1 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 0 0
    0 0 0 0 0 0 g 0 0 0 0 0
    0 0 0 0 0 0 0 1 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 1 0 0
    0 0 0 0 0 0 0 0 0 0 -g 0
    0 0 0 0 0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 0 0 0 0 0];

Bbar0 = [
    0 0 0 0
    1/mass 0 0 0
    0 0 0 0
    0 1/Ixx 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 1/Iyy 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 1/Izz];

perm0 = [2 4 8 12 6 10 1 3 7 11 5 9];
Ahatc = Ahat0(perm0, perm0);
Bbarc = Bbar0(perm0,:);
C = zeros(1,n); C(1) = 1;
Ts = 0.1; % time scale of the discritization
sysd = c2d(ss(Ahatc, Bbarc, C, 0), Ts);

Ahat = sysd.A; Ahat(abs(Ahat)<1e-6) = 0;
Bbar = sysd.B; Bbar(abs(Bbar)<1e-6) = 0;

Bhat = [Bbar, zeros(n,n-m)];
Qhat = -1*ones(Nc) + Nc*eye(Nc);

Lhat = [kron(eye(4), [1; -1]), zeros(8)];
% matrices to decouple
A = kron(eye(Nc), Ahat);
B = kron(eye(Nc), Bhat);
R = eye(n*Nc);
Q = kron(Qhat, eye(n));

DG = diag(rand(n*Nc*2, 1));
G = DG*kron(eye(n*Nc), [1; -1]);

DL = diag(rand(m*Nc*2, 1));
L = DL*kron(eye(Nc), Lhat);

%% SBD
MatA = {A, B, Q, A', B'};
MatG = {G, L};

NN = 10;
tol = 1e-6;
delta = 1e-6;
[T,X, info] = GSBD_new(MatA,MatG,NN,tol,delta);

% [T, Xm, ~] = commdec(MatA, 0);
%% applying transformation

At = T'*A*T; At(abs(At)<delta) = 0;
Bt = T'*B*T; Bt(abs(Bt)<delta) = 0;
Qt = T'*Q*T; Qt(abs(Qt)<delta/100) = 0; Qt = (Qt + Qt')/2;
Rt = T'*R*T; Rt(abs(Rt)<delta/100) = 0; Rt = (Rt + Rt')/2;
Gt = G*T; Gt(abs(Gt)<delta) = 0;
Lt = L*T; Lt(abs(Lt)<delta) = 0;

invDG = diag(1./(diag(DG)));
invDL = diag(1./(diag(DL)));

G = invDG*G;
L = invDL*L;

Gt2 = invDG(info.mat_perm{1},info.mat_perm{1})*Gt(info.mat_perm{1},:);
Lt2 = invDL(info.mat_perm{2},info.mat_perm{2})*Lt(info.mat_perm{2},:);

Num_blocks = info.Num_blocks;
dim_blocks = info.dim_blocks;
blk_indx = info.blk_indx;
G_sizes = info.G_sizes;
G_indices = info.G_indices;
%%
figure

subplot(332),spy(At),title('$A_t$','interpreter','latex'); xlabel([]),xticks([]),yticks([])
subplot(333),spy(Bt),title('$B_t$','interpreter','latex'); xlabel([]),xticks([]),yticks([])
subplot(335),spy(Qt),title('$Q_t$','interpreter','latex'); xlabel([]),xticks([]),yticks([])
subplot(336),spy(Rt),title('$R_t$','interpreter','latex'); xlabel([]),xticks([]),yticks([])

subplot(3,3,[1 4 7]),spy(Gt2),title('$G_t$','interpreter','latex'); xlabel([])
subplot(3,3,[8 9]),spy(Lt2),title('$L_t$','interpreter','latex'); xlabel([])
