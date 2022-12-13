% Simultaneous block-diagonalization algorithm
% based on Maehara and Murota's method
%
% Code by F. Vides & Amirhossein Nazerian
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
% G = blkdiag(X1,X2,X3,X4,X5);
% A = kron(diag(randn(1,2)),A);
% B = kron(diag(randn(1,2)),B);
% G = kron(diag(randn(1,2)),G);
% W = orth(randn(length(A)));
% A = W*A*W.';
% B = W*B*W.';
% G = G*W.';
% [T,X, info] = GSBD({A,B},{G},10,5e-5,5e-4);
% K = abs(T.'*A*T)>1e-4;
% subplot(131),spy(K.*(T.'*A*T))
% K = abs(T.'*B*T)>1e-4;
% subplot(132),spy(K.*(T.'*B*T))
% Gt = G*T;
% Gt = Gt(info.mat_perm{1},:);
% K = abs(Gt)>1e-4;
% subplot(133),spy(K.*Gt)
%
function [T,X, info] = GSBD(A,G,N,tol,delta,sm)
    if nargin < 6
        sm = 0;
    end
    n = size(A{1},1);
    la = length(A);
    lg = length(G);
    E = speye(n);
    K0 = A{1};
    Sum = rand*A{1};
    K0 = kron(E,K0)-kron(K0.',E);
    K = K0.'*K0;
    K0 = A{1}.';
    K0 = kron(E,K0)-kron(K0.',E);
    K = K + K0.'*K0;
    for k = 2:la
        K0 = A{k};
        Sum = Sum + rand*K0;
        K0 = kron(E,K0)-kron(K0.',E);
        K = K + K0.'*K0;
        K0 = A{k}.';
        K0 = kron(E,K0)-kron(K0.',E);
        K = K + K0.'*K0;
    end
    for k = 1:lg
        K0 = G{k}.'*(G{k});
        Sum = Sum + rand*K0;
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

    Lq = zeros(n);
    for i = 1:la+lg
        if i <= la
            At = T.'*A{i}*T; At(abs(At)<delta) = 0;
            Lq = ~~(Lq + ~~At);
        else
            Gt = G{i-la}*T; Gt(abs(Gt)<delta) = 0;
    
            for ii = 1:size(Gt,1)
                [~, col] = find(Gt(ii,:));
                for j = 1:length(col)
                    for k = j+1:length(col)
                        Lq(col(j), col(k)) = 1;
                        Lq(col(k), col(j)) = 1;
                    end
                end
            end
    
        end
    end
    Adj = double(~~Lq);
    [bins,~] = conncomp(graph(Adj));
    [~, perm_T] = sort(bins);
    
    T = T(:,perm_T);
    
    [info.Num_blocks,info.dim_blocks, info.blk_indx, info.mat_perm, info.G_sizes, info.G_indices] = Dim_blocks(Adj(perm_T, perm_T), G, T, delta);

end
%==========================================================
% Complementary block-diagonalization of one sided
% transformed matrices via permutations' identification
%==========================================================
function [Num_blocks,dim_blocks, blk_indx, mat_perm, G_sizes, G_indices] = Dim_blocks(Adj, MatG, T, tol)

    [~,dim_blocks] = conncomp(graph(Adj));
    Num_blocks = length(dim_blocks);
    
    if any(dim_blocks<0)
        error('dim_blocks is wrong!')
    end
    blk_indx = nan(Num_blocks, 2);
    mat_perm = cell(1, length(MatG));
    for c=1:Num_blocks
        d = sum(dim_blocks(1:c));
        blk_indx(c,:) = [d-dim_blocks(c)+1, d];
    end
    G_indices = cell(size(MatG));
    G_sizes = G_indices;
    for i = 1:length(MatG)
        p = [];
        matrix = MatG{i}*T; matrix(abs(matrix)<tol)=0;
        matrix = logical(matrix);
        G_indices{i} = nan(Num_blocks, 2);
        G_sizes{i} = nan(Num_blocks, 2);
        for c = 1:Num_blocks
            mat = matrix(:,blk_indx(c,1):blk_indx(c,2));
            mat = sum(mat,2);
            [I, ~] = find(mat);
            p = [p;I]; %#ok
    
            G = matrix(:,blk_indx(c,1):blk_indx(c,2));
            g = sum(G, 2);
            [~, I] = find(g);
            G_sizes{i}(c,:) = [length(I), dim_blocks(c)];
            d = sum(G_sizes{i}(1:c, 1));
            G_indices{i}(c,:) = [d-G_sizes{i}(c,1)+1, d];
        end
        p = unique(p, 'stable');
        rg = size(matrix, 1);
        if length(p) == rg
            mat_perm{i} = p;
        else
            C = setdiff(1:rg,p);
            p = [p;C']; %#ok
            if isempty(setdiff(1:rg,p)) && (length(p)==rg)
                mat_perm{i} = p;
            else
                error('Permutation is wrong!')
            end
        end
    
    end

end

