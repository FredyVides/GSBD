function [Adj, nw] = HVAC_Adj(nz, C, Cf, Ca)

snz = sqrt(nz);

nw = 4*nz + 2*snz;
wf = 2*nz;
wa = 4*snz;
% wz = 2*(nz - snz);

Adj = zeros(nz+nw+1);
% Nodes order:
% 1 to wf: floor/ceiling walls
% wf+1 to wf+wa: ambient walls
% wf+wa+1 to nw: internal walls
% nw+1 to nw+nz: zones
% nw+nz+1: ambient

% connecting floor/ceiling nodes to the zones
for i = 1:nz
    Adj(2*i-1:2*i,nw+i) = Cf*ones(2,1);
end
% connecting floor/ceiling nodes to the ambient
Adj(1:wf,end) = Cf*ones(wf,1);

% connecting ambient walls to the zones
grid = reshape(nw+1:nw+nz, [snz, snz])';
grid1 = grid; % saving it for later
if snz > 2
    for i = 2:snz-1
        for j = i:snz-1
            grid(i,j) = 0;
            grid(j,i) = 0;
        end
    end
end
gridv = reshape(grid', [length(grid)^2, 1]);
gridv(gridv == 0) = [];
corners = [min(gridv), min(gridv)+snz-1, max(gridv)-snz+1, max(gridv)];
j = 1;
for i = 1:length(gridv)
    Adj(wf+j,gridv(i)) = Ca;
    j = j + 1;
    if any(gridv(i) == corners)
        Adj(wf+j,gridv(i)) = Ca;
        j = j + 1;
    end
end

% connecting ambient walls to the ambient
Adj(wf+1:wf+wa,end) = Ca*ones(wa,1);

% connecting internal walls to the zones
k = 1;
for i = 1:snz
    for j = 1:snz
        if i+1 <= snz
            Adj(wf+wa+k, grid1(i,j)) = C;
            Adj(wf+wa+k, grid1(i+1,j)) = C;
            k = k + 1;
        end
        if j+1<=snz
            Adj(wf+wa+k, grid1(i,j)) = C;
            Adj(wf+wa+k, grid1(i,j+1)) = C;
            k = k + 1;
        end
    end
end
Adj = Adj + Adj';

end