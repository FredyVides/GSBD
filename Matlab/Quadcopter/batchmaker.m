function [Abar, Bbar, Qbar, Rbar] = batchmaker(A, B, Q, R, Nh, n)

Bbar = zeros(n, Nh*n);
for i = 0:Nh-1
    temp = [];
    counter = i:-1:0;
    for j = counter
        temp1 = A^j*B;
        temp = [temp, temp1];
    end
    temp = [temp, zeros(n, (Nh-length(counter))*n)];
    Bbar = [Bbar; temp];
end
Rbar = kron(eye(Nh), R);
Qbar = kron(eye(Nh+1), Q);
Abar = eye(n);
for i = 1:Nh
    Abar = [Abar;A^i];
end