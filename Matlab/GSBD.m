% Generalized Simultaneous block-diagonalization (GSBD) algorithm
%
% Code by F. Vides &
% A. Nazerian
% 
% Date of current version: March 25, 2023
%
% Ref: 
% A. Nazerian, F. Vides and F. Sorrentino, "Algebraic Decomposition 
% of Model Predictive Control Problems," in IEEE Control Systems Letters, 
% vol. 7, pp. 1441-1446, 2023, doi: 10.1109/LCSYS.2023.3252162.
%
%% inputs: 
% A: cell of square matrices for transformation of both sides:
% At{i} = T'*A{i}*T;
%
% G: cell of rectangular matrices for transformation of one sides: 
% Gt{i} = G{i}*T; Gt{i} = Gt{i}(info.mat_perm{i}, :);
%
% N: number of eigenvalues for 'eigs'
% tol: tolerance for 'eigs'
% delta: tolerance to set the eintries of the transformed matrix to zero 
% sprs: run the code using the sparse matrix mode or the full matrix mode
%       for full matrices, sprs=0
%       for sparse matrices, sprs=1
%
%% outputs:
% T: square transformation matrix
%
% X: commutant matrix
%
% info: struct with information about:
% 1: number of blocks after transformation
% 2: dimensions of the blocks
% 3: indices of the blocks for set of transformed At{:} matrices
% 4: permutation for set of transformed Gt{:} matrices
% 5: sizes of the blocks of the permuted Gt{:}
% 6: indexes of the blocks (rows) of the permuted Gt{:}

%% 
function [T,X, info] = GSBD(A,G,N,tol,delta, sprs)
    
    n = size(A{1},1);
    la = length(A);
    lg = length(G);
    
    if sprs
        for i = 1:la+lg
            if i<=la
                A{i} = sparse(A{i});
            else
                G{i-la} = G{i-la};
            end
        end
    end
    
    E = speye(n);
    K0 = A{1};
    K0 = kron(E,K0)-kron(K0.',E);
    K = rand*(K0.'*K0);
    K0 = A{1}.';
    K0 = kron(E,K0)-kron(K0.',E);
    K = K + rand*(K0.'*K0);
    for k = 2:la
        K0 = A{k};
        K0 = kron(E,K0)-kron(K0.',E);
        K = K + rand*(K0.'*K0);
        K0 = A{k}.';
        K0 = kron(E,K0)-kron(K0.',E);
        K = K + rand*(K0.'*K0);
    end
    for k = 1:lg
        if sprs
            D = diag(sparse(rand(size(G{k}, 1), 1)));
        else
            D = diag(rand(size(G{k}, 1), 1));
        end
        K0 = G{k}.'*D*(G{k});
        K0 = kron(E,K0)-kron(K0.',E);
        K = K + rand*(K0.'*K0);
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
    X(abs(X)<delta) = 0;
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
    Adj = double(~~Lq+Lq');
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

