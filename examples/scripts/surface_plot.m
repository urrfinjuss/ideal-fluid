N = 4000;
u = pi*(2*(0:N-1)/N - 1.);
a1 =  0.25*pi + 0.1i;
b1 =  0.25*pi + 0.2i;

a2 = -conj(a1);  
b2 = -conj(b1);

%-------------
A = -2.2385i;

dip1 = (log(i*(sin(0.5*(u - a1)))) - log(i*sin(0.5*(u - b1))));
dip2 = (log(i*(sin(0.5*(u - a2)))) - log(i*sin(0.5*(u - b2))));
z = u + A*dip1 + A*dip2;

fh = fopen('splash_001.txt','w');
fprintf(fh, '# 1. u 2. x 3. y\n\n');
for j = 1:N
  fprintf(fh, '%.12e\t%.12e\t%.12e\n', u(j), real(z(j)), imag(z(j)));
end
fclose(fh);
