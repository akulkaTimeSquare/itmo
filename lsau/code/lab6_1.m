n = 6;
k = 5;
tp = 8;
alphas = [1, ];
for j = 1:n
    p = [1, 1];
    alphas = conv(alphas, p);
end
alphas

a = zeros(n);
tpn = 12;
omega = tpn/tp;
for j = 1:n
    a(j) = alphas(j)*omega^(n-j+1);
end
a = cat(1, [1, ], wrev(a(:, 1)))
b = a(end)*k
r = roots(a)
rp = real(r);
ip = imag(r);
nu = abs(rp(3));
t = 1/nu*log(1/0.05)