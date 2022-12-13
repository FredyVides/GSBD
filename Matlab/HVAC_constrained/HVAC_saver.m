clc
clear
close all

%% matrices from quadcopter
for zone = (2:6).^2 % number of quadcopters (variable)
    clearvars -except zone
    close all
    %     data = load(['zone', num2str(zone)]);
    %     n = data.n; % number of states
    %     m = data.m; % number of inputs
    %
    %     % matrices to decouple
    %     A = data.Ad; A(abs(A)<1e-6) = 0;
    %     B = data.Bd; B(abs(B)<1e-6) = 0;
    %     R = data.R;
    %     Q = data.Q;
    [A, B, Q, R] = HVACcreator(zone, 1e-6, 600);
    n = length(A);
    m = zone;
    %     DG = diag(rand(length(data.Aww)*2, 1));
    %     G = DG*kron(eye(length(data.Aww)), [1; -1]);
    %     G = [G, zeros(size(G,1), n-size(G,2))];

    G = kron(eye(zone), [1; -1]);
    DG = diag(rand(size(G,1), 1));
    G = DG*[zeros(size(G,1), n-size(G,2)), G];

    %     G = kron(eye(n), [1; -1]);
    %     DG = diag(rand(size(G,1), 1));
    %     G = DG*G;

    DL = diag(rand(m*2, 1));
    Lhat = [kron(eye(m), [1; -1]), zeros(2*m, n-m)];
    L = DL*Lhat;

    %% SBD
    MatA = {sparse(A), sparse(B), sparse(Q), sparse(A'), sparse(B'), sparse(R)};
    MatG = {G, L};

    NN = 20;
    tol = 1e-8;
    delta = 1e-6;
    tic
    [T,X, info] = GSBD_new(MatA,MatG,NN,tol,delta);
    toc

    % [T, Xm, ~] = commdec(MatA, 0);
    %% applying transformation

    At = T'*A*T; At(abs(At)<delta) = 0;
    Bt = T'*B*T; Bt(abs(Bt)<delta) = 0;
    Qt = T'*Q*T; Qt(abs(Qt)<delta) = 0; Qt = (Qt + Qt')/2;
    Rt = T'*R*T; Rt(abs(Rt)<delta) = 0; Rt = (Rt + Rt')/2;
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
    dim_blocks
    %%
    if issymmetric(Qt) && issymmetric(Rt) && size(G, 1) == size(Gt2, 1) && ...
            size(L, 1) == size(Lt2, 1) && sum(dim_blocks) == n && ...
            norm(sort(eig(A))-sort(eig(At))) < 10*delta
        disp('save it')
%         save(['newzone', num2str(zone)])
    end

end