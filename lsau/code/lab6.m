n = 6;
k = 5;
tp = 8;
alphas = [1, ];
for j = 1:n
    p = [1, -exp(1i*(pi/2 + (2*j-1)/(2*n)*pi))];
    alphas = conv(alphas, p);
end
alphas = real(alphas)

a = zeros(n);
tpn = 14.1;
omega = tpn/tp;
for j = 1:n
    a(j) = alphas(j)*omega^(n-j+1);
end
a = cat(1, [1, ], wrev(a(:, 1)))
b = a(end)*k
r = roots(a);
rp = real(r)
ip = imag(r)
nu = abs(rp(1));
t = 1/nu*log(1/0.05)
