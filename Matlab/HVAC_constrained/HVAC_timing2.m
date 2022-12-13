clc
clear
close all
%% variables
zones = (2:6).^2;
Horizons = [10 20 30];
%% results
TimeOrg = nan(length(zones), length(Horizons));
STD_org = TimeOrg;
TimeTran = TimeOrg;
STD_trans = TimeOrg;
%% first loops
for izone = 1:length(zones)
    zone = zones(izone); % number of quadcopters (variable)

    clearvars -except zones Horizons TimeOrg STD_org TimeTran STD_trans zone Nh izone
    % loading data
    data = load(['newzone', num2str(zone)]);

    % original domain
    org.A = data.A;
    org.B = data.B;
    org.Q = data.Q;
    org.R = data.R;
    org.G = data.G;
    org.L = data.L;
    g = ones(size(org.G, 1), 1);
    l = ones(size(org.L, 1), 1);

    % transfomed domain
    tran.A = data.At;
    tran.B = data.Bt;
    tran.Q = data.Qt;
    tran.R = data.Rt;
    tran.GG = data.Gt2;
    tran.LL = data.Lt2;
    gt = g(data.info.mat_perm{1});
    lt = l(data.info.mat_perm{2});


    n = data.n;
    blk_indx = data.blk_indx;
    dim_blocks = data.dim_blocks;

    twoblock = true; %%%%%%%%%% attaches all scalar blocks to eachother
    if twoblock
        [~, I1] = find(dim_blocks==1);
        [~, IL] = find(dim_blocks>1);
        perm = [];
        for i = [IL, I1]
            perm = [perm, blk_indx(i, 1):blk_indx(i,2)];
        end
        tran.A = tran.A(perm,perm);
        tran.B = tran.B(perm,perm);
        tran.Q = tran.Q(perm,perm);
        tran.R = tran.R(perm,perm);
        tran.GG = tran.GG(:,perm);
        tran.LL = tran.LL(:,perm);
        Num_blocks = length(IL) + 1;
        dim_blocks = [dim_blocks(IL), length(I1)];
        blk_indx = nan(Num_blocks, 2);
        for c=1:Num_blocks
            d = sum(dim_blocks(1:c));
            blk_indx(c,:) = [d-dim_blocks(c)+1, d];
        end
        [G_sizes, G_indices] = Gsize(tran.GG, Num_blocks, dim_blocks, blk_indx);
        [L_sizes, L_indices] = Gsize(tran.LL, Num_blocks, dim_blocks, blk_indx);
    else
        Num_blocks = data.Num_blocks;

        [G_sizes, G_indices] = Gsize(tran.GG, Num_blocks, dim_blocks, blk_indx);
        [L_sizes, L_indices] = Gsize(tran.LL, Num_blocks, dim_blocks, blk_indx);
    end
    %% Second loop
    for iNh = 1:length(Horizons) % Horizon (variable)
        clc
        disp(['zone = ', num2str(zone)])
        Nh = Horizons(iNh);
        % Batch Implementation
        disp(['Nh = ', num2str(Nh)])
        % original domain
        [org.Ab, org.Bb, org.Qb, org.Rb] = batchmaker(org.A, org.B, org.Q, org.R, Nh, n);

        % transformed domain
        tran.Ab = cell(Num_blocks,1);
        tran.Bb = tran.Ab; tran.Qb = tran.Ab; tran.Rb = tran.Ab;  tran.g = tran.Ab;
        tran.l = tran.Ab; tran.G = tran.Ab; tran.L = tran.Ab;

        for i = 1:Num_blocks
            [tran.Ab{i}, tran.Bb{i}, tran.Qb{i}, tran.Rb{i}] = ...
                batchmaker(tran.A(blk_indx(i,1):blk_indx(i,2), blk_indx(i,1):blk_indx(i,2)), ...
                tran.B(blk_indx(i,1):blk_indx(i,2), blk_indx(i,1):blk_indx(i,2)), ...
                tran.Q(blk_indx(i,1):blk_indx(i,2), blk_indx(i,1):blk_indx(i,2)), ...
                tran.R(blk_indx(i,1):blk_indx(i,2), blk_indx(i,1):blk_indx(i,2)), ...
                Nh, dim_blocks(i) ...
                );
            tran.g{i} = gt(G_indices(i,1):G_indices(i,2));
            tran.G{i} = tran.GG(G_indices(i,1):G_indices(i,2), blk_indx(i,1):blk_indx(i,2));

            tran.l{i} = lt(L_indices(i,1):L_indices(i,2));
            tran.L{i} = tran.LL(L_indices(i,1):L_indices(i,2), blk_indx(i,1):blk_indx(i,2));
        end

        %% solving constrained MPC
        maxiter = 1; % parameter
        options = optimoptions('quadprog','Display','none', 'Algorithm','interior-point-convex', 'MaxIterations', 5000);

        % original
        H = 2*(org.Rb + org.Bb'*org.Qb*org.Bb);
        H(abs(H)<1e-8) = 0;
        H = (H'+H)/2;

        Hx = kron(eye(Nh+1), org.G);
        Kx = kron(ones(Nh+1,1), g);
        Hu = kron(eye(Nh), org.L);
        Ku = kron(ones(Nh,1), l);
        A = [Hx*org.Bb; Hu];
        A(abs(A)<1e-8) = 0;
        % H = sparse(H); A = sparse(A);

        if issparse(H)
            mode = true;
        else
            mode = false;
        end

        X0 = 0.1*(2*rand(size(org.Ab, 2), maxiter) - 1);
        Z0 = nan(size(X0));
        t_org = nan(maxiter, 1);
        for iter = 1:maxiter
            x0 = X0(:,iter);
            Z0(:,iter) = ((data.T)')*X0(:,iter);
            f = 2*org.Bb'*org.Qb'*org.Ab*x0;
            b = [Kx - Hx*org.Ab*x0; Ku];

%             exitflag = 0;
%             while ~(exitflag == 1)
                % tStart = cputime;
                % tic
                % [U_org,~,exitflag,~] = quadprog(H,f,A,b,[],[],[],[], [], options);
                % tEnd = cputime - tStart;
                % tEnd = toc;
                fun = @() quadprog(H,f,A,b,[],[],[],[], [], options);
                tEnd = timeit(fun);
%                 if ~(exitflag == 1)
%                     x0 = 0.2*(2*rand(size(org.Ab, 2), 1) - 1);
%                     X0(:,iter) = x0;
%                     Z0(:,iter) = data.T'*X0(:,iter);
%                     f = 2*org.Bb'*org.Qb'*org.Ab*x0;
%                     b = [Kx - Hx*org.Ab*x0; Ku];
%                 end
%             end
            t_org(iter) = tEnd;

        end

        TimeOrg(izone, iNh) = mean(t_org);
        STD_org(izone, iNh) = std(t_org);
        %         Sol_org = reshape(U_org, n, Nh);
        %% transformed
        t_trans = nan(maxiter, Num_blocks);

        %         U_tran = [];
        if  mode
            % solving in parallel
            for i = 1:Num_blocks
                H = round(2*(tran.Rb{i} + ((tran.Bb{i})')*tran.Qb{i}*tran.Bb{i}), 8);
                H = (H'+H)/2;
                Hx = kron(eye(Nh+1), tran.G{i});
                Kx = kron(ones(Nh+1,1), tran.g{i});
                Hu = kron(eye(Nh), tran.L{i});
                Ku = kron(ones(Nh,1), tran.l{i});
                A = [Hx*tran.Bb{i}; Hu];

                time_trans = nan(maxiter, 1);
                for iter = 1:maxiter
                    z0 = Z0(blk_indx(i,1):blk_indx(i,2),iter);
                    if any(tran.G{i}*z0 - tran.g{i} > 0)
                        error('initial condition is infeasible')
                    end
                    f = 2*tran.Bb{i}'*tran.Qb{i}'*tran.Ab{i}*z0;
                    b = [Kx - Hx*tran.Ab{i}*z0; Ku];

                    %tStart = cputime;
                    tic
                    [Ui_tran,~,exitflag,~] = quadprog(H,f,A,b,[],[],[],[], [], options);
                    %tEnd = cputime - tStart;
                    tEnd = toc;
                    if exitflag == 1
                        time_trans(iter) = tEnd;
                        %                     U_tran = [U_tran; Ui_tran];
                    else
                        error(['Exitflag of trans! = ', num2str(exitflag)])
                    end
                end
                t_trans(:, i) = time_trans;
            end
            sum_Ttrans = max(t_trans, [], 2);
        else
            % solving in series
            for i = 1:Num_blocks
                H = round(2*(tran.Rb{i} + ((tran.Bb{i})')*tran.Qb{i}*tran.Bb{i}), 8);
                H = (H'+H)/2;
                Hx = kron(eye(Nh+1), tran.G{i});
                Kx = kron(ones(Nh+1,1), tran.g{i});
                Hu = kron(eye(Nh), tran.L{i});
                Ku = kron(ones(Nh,1), tran.l{i});
                A = [Hx*tran.Bb{i}; Hu];
                
                if i == 2
                   H = sparse(H); A = sparse(A);
                end

                time_trans = nan(maxiter, 1);
                for iter = 1:maxiter
                    z0 = Z0(blk_indx(i,1):blk_indx(i,2),iter);
                    if any(tran.G{i}*z0 - tran.g{i} > 0)
                        error('initial condition is infeasible')
                    end
                    f = 2*tran.Bb{i}'*tran.Qb{i}'*tran.Ab{i}*z0;
                    b = [Kx - Hx*tran.Ab{i}*z0; Ku];

                    % tStart = cputime;
                    % tic
                    % [Ui_tran,~,exitflag,~] = quadprog(H,f,A,b,[],[],[],[], [], options);
                    % tEnd = cputime - tStart;
                    % tEnd = toc;
                    fun = @() quadprog(H,f,A,b,[],[],[],[], [], options);
                    tEnd = timeit(fun);
%                     if exitflag == 1
                        time_trans(iter) = tEnd;
%                         %                     U_tran = [U_tran; Ui_tran];
%                     else
%                         error(['Exitflag of trans! = ', num2str(exitflag)])
%                     end
                end
                t_trans(:, i) = time_trans;
            end
            sum_Ttrans = sum(t_trans, 2);
        end

        TimeTran(izone, iNh) = mean(sum_Ttrans);
        STD_trans(izone, iNh) = std(sum_Ttrans);

        %         Sol_tran = reshape(U_tran, n, Nh);
        %         norm(data.T*Sol_tran - Sol_org)
    end
end

if mode
    modetxt = 'sparse';
else
    modetxt = 'full';
end
txt = ['HVAC_zone', num2str(min(zones)),'_',num2str(max(zones)),...
    '_H',num2str(min(Horizons)),'_',num2str(max(Horizons)),'_', modetxt, 'blocks2.mat'];

save(txt)
%% plots

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
hAx=axes;
hAx.YScale='log';
hold all
shape = {'d', 'o', 'p', 'square'};
blue = [0.00,0.45,0.74];
red = [0.85,0.33,0.10];
for i = 1:length(Horizons)
    if i == 2 && length(Horizons) == 4
        continue
    end
    plot(sqrt(zones), TimeOrg(:,i), Color=blue,  LineStyle=':', Marker = shape{i},  LineWidth=3, MarkerSize=15, MarkerFaceColor='auto', DisplayName=['$N = $ ', num2str(Horizons(i)), ', Org.'])
    hold on
    plot(sqrt(zones), TimeTran(:,i), Color=red,  LineStyle=':', Marker = shape{i},  LineWidth=3, MarkerSize=15, MarkerFaceColor='auto', DisplayName=['$N = $ ', num2str(Horizons(i)), ', Trans.'])
end
set(gca, 'fontname', 'Times New Roman', 'fontsize', 40)
xlabel('$N_c$', 'Interpreter','latex')
ylabel('CPU time (seconds)', 'Interpreter','latex')
% legend('Original', 'Transformed', 'interpreter', 'latex', 'location', 'Orientation','horizontal')
lgnd = legend('show');
lgnd.Interpreter = 'latex';
lgnd.Location = 'northwest';
lgnd.NumColumns=4;
lgnd.FontSize=30;
xticks(sqrt(zones))
% xticklabels({'$2^2$', '$3^2$', '$4^2$', '$5^2$', '$6^2$'})
xlim([1, max(sqrt(zones))+1])
yticks([1e-0 1e1 1e2 1e3 1e4])
% yticklabels({'0.001', '0.01', '0.1', '1'})
% ylim([0.001 10])
set(gca,'TickLength',[0.03,0.025], 'layer','top', 'XGrid', 0, 'YGrid', 0, 'LineWidth', 2, 'Box','on');
%%

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
hAx=axes;
hAx.YScale='log';
hold all
shape = {'d', 'o', 'p', 'square'};
blue = [0.00,0.45,0.74];
red = [0.85,0.33,0.10];
for i = 1:length(Horizons)
    errorbar(sqrt(zones), TimeOrg(:,i), STD_org(:,i), Color=blue,  LineStyle=':', Marker = shape{i},  LineWidth=3, MarkerSize=10, MarkerFaceColor='auto', DisplayName=['$N = $ ', num2str(Horizons(i)), ', Original'])
    hold on
    errorbar(sqrt(zones), TimeTran(:,i), STD_trans(:,i), Color=red,  LineStyle=':', Marker = shape{i},  LineWidth=3, MarkerSize=10, MarkerFaceColor='auto', DisplayName=['$N = $ ', num2str(Horizons(i)), ', Transformed'])
end
set(gca, 'fontname', 'Times New Roman', 'fontsize', 30)
xlabel('Number of quadcopters', 'Interpreter','latex')
ylabel('CPU time (seconds)', 'Interpreter','latex')
% legend('Original', 'Transformed', 'interpreter', 'latex', 'location', 'Orientation','horizontal')
lgnd = legend('show');
lgnd.Interpreter = 'latex';
lgnd.Location = 'northoutside';
lgnd.NumColumns=4;
lgnd.FontSize=20;

xticks(sqrt(zones))
% xticklabels({'$2^2$', '$3^2$', '$4^2$', '$5^2$', '$6^2$'})
xlim([1, max(sqrt(zones))+1])
% yticks([1e-3 1e-2 1e-1 1])
% yticklabels({'0.001', '0.01', '0.1', '1'})
% ylim([0.001 10])
set(gca,'TickLength',[0.03,0.025], 'layer','top', 'XGrid', 0, 'YGrid', 0, 'LineWidth', 2, 'Box','on');
%%
function [G_sizes, G_indices] = Gsize(Gt, Num_blocks, dim_blocks, blk_indx)
Gt = ~~Gt;
G_indices = nan(Num_blocks, 2);
G_sizes = nan(Num_blocks, 2);
for i = 1:Num_blocks
    G = Gt(:,blk_indx(i,1):blk_indx(i,2));
    g = sum(G, 2);
    [~, I] = find(g);
    G_sizes(i,:) = [length(I), dim_blocks(i)];
    d = sum(G_sizes(1:i, 1));
    G_indices(i,:) = [d-G_sizes(i,1)+1, d];
end

end