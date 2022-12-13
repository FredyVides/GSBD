% Simultaneous block-diagonalization algorithm
% based on Maehara and Murota's method
%
% Code by F. Vides
%
% Example:
%
% X1 = randn(4);
% X2 = randn(7);
% X3 = randn(3);
% X4 = randn(5);
% X5 = randn(2);
% A = blkdiag(X1,X2,X3,X4,X5);
% X1 = randn(4);
% X2 = randn(7);
% X3 = randn(3);
% X4 = randn(5);
% X5 = randn(2);
% B = blkdiag(X1,X2,X3,X4,X5);
% X1 = randn(4);
% X2 = randn(7);
% X3 = randn(3);
% X4 = randn(5);
% X5 = randn(2);
% C = blkdiag(X1,X2,X3,X4,X5);
% A = kron(diag(randn(1,2)),A);
% B = kron(diag(randn(1,2)),B);
% C = kron(diag(randn(1,2)),C);
% W = orth(randn(length(A)));
% A = W*A*W.';
% B = W*B*W.';
% C = W*C*W.';
% [T,X] = SBD({A,B,C},10,5e-5);
% K = abs(T.'*A*T)>1e-4;
% subplot(131),spy(K.*(T.'*A*T))
% K = abs(T.'*B*T)>1e-4;
% subplot(132),spy(K.*(T.'*B*T))
% K = abs(T.'*C*T)>1e-4;
% subplot(133),spy(K.*(T.'*C*T))
%
function [T,X] = SBD(A,N,tol)
  n = size(A{1},1);
	la = length(A);
	E = speye(n);
	K0 = A{1};
	K0 = kron(E,K0)-kron(K0.',E);
	K = K0.'*K0;
        K0 = A{1}.';
	K0 = kron(E,K0)-kron(K0.',E);
	K = K + K0.'*K0;
	for k = 2:la
		K0 = A{k};
		K0 = kron(E,K0)-kron(K0.',E);
		K = K + K0.'*K0;
        	K0 = A{k}.';
		K0 = kron(E,K0)-kron(K0.',E);
		K = K + K0.'*K0;
	end
	[u,l] = eigs(K,N,tol);
	l = abs(diag(l));
	f = find((l-tol)<tol);
	lf = length(f);
	c = randn(lf,1);
	c = c/norm(c);
	X = u(:,f)*c;
	X = reshape(X,n,n);
	X = (X+X.')/2;
	[T,~] = eig(X);
end
