A = [1 0; 0 0];
b = [2; 4];
d = [2; 4];

A0f = [0 1; -0.25 -1];
b0f = [0; 1];
N = b0f/d;

t_n = 0.9;
t_n_star = 4.7;
w0 = t_n_star / t_n;
A_zh = [0 1; -w0^2 -2*w0];

H = [1 0];
M = lyap(A, -A_zh, -b*H);

K = H * inv(M);
Am = A - b*K;

Q = eye(2);
P = lyap(Am', Q);
disp(P);