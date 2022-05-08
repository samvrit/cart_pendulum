Ts = 1e-4;
dsys = c2d(sys, Ts);
F = dsys.A;
H = dsys.C;
Q = (dsys.B * dsys.B') .* 1000;
% F = A;
% H = C;
% Q = (B * B') * (Ts*Ts) * 1000;
% Q = 1000;

R = eye(4) .* 1;
P = Q;
n = 20000;
K_Kalman = zeros(4,4);
I = eye(4);

[~, K_Kalman_matlab] = kalman(dsys, 1000, 1);

for i = 1:n
    P = (F * P * F') + Q;
    S = (H * P * H') + R;
    K_Kalman = (P * H') / S;
    P = (I - (K_Kalman * H)) * P;
end
disp(K_Kalman)