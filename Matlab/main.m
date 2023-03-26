clc
clear
close all
%% Example:

% creating random block diagonal matrices
X1 = randn(4);
X2 = randn(7);
X3 = randn(3);
X4 = randn(5);
X5 = randn(2);
A = blkdiag(X1,X2,X3,X4,X5);

X1 = randn(4);
X2 = randn(7);
X3 = randn(3);
X4 = randn(5);
X5 = randn(2);
B = blkdiag(X1,X2,X3,X4,X5);

X1 = rand(6, size(X1,2));
X2 = rand(3, size(X2,2));
X3 = rand(4, size(X3,2));
X4 = rand(8, size(X4,2));
X5 = rand(2, size(X5,2));
G = blkdiag(X1,X2,X3,X4,X5);

X1 = rand(2, size(X1,2));
X2 = rand(2, size(X2,2));
X3 = rand(4, size(X3,2));
X4 = rand(5, size(X4,2));
X5 = rand(3, size(X5,2));
L = blkdiag(X1,X2,X3,X4,X5);

A = kron(diag(randn(1,2)),A);
B = kron(diag(randn(1,2)),B);
G = kron(diag(randn(1,2)),G);
L = kron(diag(randn(1,2)),L);

figure('Name', 'Original')
subplot(221),spy(A)
subplot(222),spy(B)
subplot(223),spy(G)
subplot(224),spy(L)

% perturbing the matrices by a random orthogonal transformation
W = orth(randn(length(A)));
Ap = W*A*W.';
Bp = W*B*W.';
Gp = G*W.';
Lp = L*W.';

figure('Name', 'Perturbed')
subplot(221),spy(Ap)
subplot(222),spy(Bp)
subplot(223),spy(Gp)
subplot(224),spy(Lp)

% finding the GSBD transformation
[T,X, info] = GSBD({Ap,Bp},{Gp, Lp},10,1e-8,1e-6, 0);

% retreving the original structure through the GSBD transformation
figure('Name', 'Transformed')
K = abs(T.'*Ap*T)>1e-6;
subplot(221),spy(K.*(T.'*Ap*T))

K = abs(T.'*Bp*T)>1e-6;
subplot(222),spy(K.*(T.'*Bp*T))

Gt = Gp*T;
Gt = Gt(info.mat_perm{1},:); % applying permutation
K = abs(Gt)>1e-6;
subplot(223),spy(K.*Gt)

Lt = Lp*T;
Lt = Lt(info.mat_perm{2},:); % applying permutation
K = abs(Lt)>1e-6;
subplot(224),spy(K.*Lt)
