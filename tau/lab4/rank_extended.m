A = [11 -2 13; 6 -1 6; -6 -1 -8];
Bf = [-6 0 0 -1; 0 0 0 0; 6 0 0 0];
G = [35 56 22 -42; -11 -17 -7 12; -6 -10 -5 10; 11 18 6 -13];
C = [2 -2 1];
D = [1 2 1 -1];

Axw = [A Bf; zeros(4, 3) G];
Cxw = [C D];

V = [Cxw; Cxw*Axw; Cxw*Axw^2; Cxw*Axw^3; Cxw*Axw^4; Cxw*Axw^5; Cxw*Axw^6]
rank(V)