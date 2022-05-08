Ts = 1e-4;

t = 0:Ts:9.9999;

N = size(t, 2);

x = zeros(1, N);
u1 = zeros(1, N);
theta = zeros(1, N);
v1 = zeros(1, N);

F = zeros(1, N);

m1 = 1;  % kg
m2 = 0.5;  % kg
l = 5;  % m
b = 0.0;  % Nm/rad/s
g = 9.81;  % m/s^2

theta(1) = deg2rad(10);

A = [0 1 0 0;
     0 0 0 0;
     0 0 0 1;
     0 0 g/l 0];
B = [0; 1/m1; 0; 1/(m1*l)];

C = diag([1 0 1 0]);

D = zeros(4, 1);

Q = diag([1 10 1000 1000]);
R = 1;

K = lqr(A, B, Q, R);

sys = ss(A, B, C, D);

Q_kalman = 1000;
R_kalman = 1;

[Kf, L] = kalman(sys, Q_kalman, R_kalman, 0);

q = zeros(4, N);

x_hat = zeros(4, N);

for i = 2:N
    
    F(i) = -K * x_hat(:,i-1);
    
    c1 = (m1 + m2 - m2*cos(theta(i-1))*cos(theta(i-1)));
    term1 = F(i);
    term2 = m2 * v1(i-1) * v1(i-1) * l * sin(theta(i-1));
    term3 = m2 * g * sin(theta(i-1)) * cos(theta(i-1));
    term4 = (b / l) * v1(i-1) * cos(theta(i-1));
    u1(i) = u1(i-1) + (Ts/c1) * (term1 - term2 + term3 - term4);
    x(i) = x(i-1) + Ts*u1(i-1);
    
    term1 = cos(theta(i-1)) * F(i) / l;
    term2 = m2 * v1(i-1) * v1(i-1) * sin(theta(i-1)) * cos(theta(i-1));
    term3 = g * (m1+m2) * sin(theta(i-1)) / l;
    term4 = b * (m1+m2) * v1(i-1) / (m2 * l * l);
    v1(i) = v1(i-1) + (Ts/c1) * (term1 - term2 + term3 + term4); 
    theta(i) = theta(i-1) + Ts*v1(i-1);
    
    q(:,i) = [x(i); u1(i); theta(i); v1(i)];
    
    y = C*q(:,i);
    
    x_hat(:,i) = x_hat(:,i-1) + Ts*( (A-B*K)*x_hat(:,i-1) + L*(y - C*x_hat(:,i-1)) );
end

% plot(t, q, t, x_hat);

