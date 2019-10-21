function [out, ly, ny, M] = PetviashviliMethodMk5(Nin, Lambda, Omega_in, A, k0, nsteps)
% this code doesn't work, and I suspect this is 
% because y_u/|z_u| is computed not very accurately.
% Two things:  -- do iterations in k-space
%	       -- write z_u/|z_u| with Fourier series
% N = 16384; z = PetviashviliMethod(N, 1.06, 0.05, 4, 1000);
% N = 65536; z = PetviashviliMethod(N, 1.2, 0.25, 4, 4000);  
global k Ak Pk Hk Omega Sigma



  graphics_toolkit("gnuplot")
  k = fftshift(-Nin/2:Nin/2-1);
  u = pi*(2*(0:Nin-1)/Nin - 1);
  Ak = -1i./k;
  Ak(1) = 0;
 
  Pk = 0.5*(1-sign(k));
  Hk = 1.i*sign(k);
  
  Sigma = 1;
  Omega = Omega_in;
  zk = zeros(1,Nin);
  zk(end) = -1i;
  zk(end-k0) = 1i*A;
  
  z  = ifft(zk)*Nin;
  yk = real(fft(0.5i*(conj(z)-z))/Nin);
  mu = 0.5*sum(abs(k).*abs(zk).^2); M = 1./sqrt(mu);
  yk = yk/sqrt(mu);
  yk0 = yk;
  for j = 1:nsteps
    M = (linear(yk, Lambda)*yk')/(nonlinear(yk)*yk');
    
    yk = M*nonlinear(yk)./(abs(k)-Lambda);
    %yk = yk/abs(2*yk(end));
    %yk = yk/sqrt(mu);
    %yk = yk/sqrt(mu);
    %yk = yk*M;
    zk = -1.*(Hk - 1.i).*yk;
    mu = 0.5*sum(abs(k).*abs(zk).^2);
    %yk = yk/sqrt(mu);
    
    %zk = -1.*(Hk - 1.i).*yk;
    %mu = 0.5*sum(abs(k).*abs(zk).^2);
    fprintf('Petviashvili iteration %4d: M = %.12e\tmu = %.12e\n', j, M, mu);
  end
  res = linear(yk,Lambda)-nonlinear(yk);
  ly = linear(yk, Lambda);
  ny = nonlinear(yk);
  NormR = sqrt(sum(abs(res).^2));
  printf('Norm of Residual = %.12e\n', NormR);
  figure(1)
  plot(u, ifft(yk0)*Nin, u, ifft(yk)*Nin)
  
  zk = -1.*(Hk - 1.i).*yk;
  out = ifft(zk)*Nin;
  
  % prepare wave for the dynamic code:
  dz = ifft(1i*k.*fft(out));
  dPhi = 1i*Omega*sqrt(-M).^2*ifft(1i*k.*Pk.*fft(abs(out).^2));
  
  R = 1./dz;
  V = 1i.*R.*dPhi;
  % write data to input file:
  fh = fopen('../config/ntravel_001.txt','w');
  fprintf(fh, '1. u 2.-3. R 4.-5. V\n');
  fprintf(fh, '# Time = %.16e\tus = %.16e\tqs = %.16e\tl = %.16e\n\n', 0.0, 0.0, 0.0, 1);
  for j = 1:Nin
    q = pi*(2*(j-1)/Nin - 1);
    fprintf(fh, '%.16e\t', q);
    fprintf(fh, '%.16e\t%.16e\t', real(R(j)), imag(R(j)));
    fprintf(fh, '%.16e\t%.16e\n', real(V(j)), imag(V(j)));
  end
  fclose(fh);


  return

end

function outk = linear(ink, Lambda)
 global k
 
 outk = real((abs(k)-Lambda).*ink);
end

function outk = nonlinear(ink)
 global k Ak Pk Hk Omega Sigma

 N = length(k);
 zk = -1.*(Hk - 1.i).*ink;
 z = ifft(zk)*N;
 dz = ifft(1i*k.*zk)*N;
 tmp = 1 + 1./abs(dz);
 %Tmp = ifft(ratio_fourier( ifftshift(fft(1.+0.*k)/N), ifftshift(fft(abs(dz)))))*N;
 %tmp = 1 - Tmp;
 kAbsZ2 = ifft(abs(k).*fft(abs(z).^2));


 outk = 0.5*fft( real(dz).*tmp - ifft(Hk.*fft(imag(dz).*tmp )))/N; 
 outk = outk + 0.25*Omega^2/Sigma.*Ak.*fft( real(z).*kAbsZ2 - ifft(Hk.*fft(imag(z).*kAbsZ2)))/N;
 outk = real(outk);
 outk(N/4:3*N/4) = 0; 
 %idk = zeros(1, N); idk(1) = 1;
 %tmp = ratio_fourier( ifftshift(idk), ifftshift(fft(abs(dZ))));
 
end

function sk = ratio_fourier(fk,gk)
  % takes Fourier coefficients fk and gk 
  % and builds sk = FT(f/g) via series inversion
  global k
  N = length(k);
  A = zeros(N+1,N+1);
  Gp = zeros(N+1,1); 
  Gp(2:N+1) = gk(N:-1:1);
  Gp(1) = gk(1);

  Fp = zeros(N+1,1);
  Fp(1:N) = fk(1:N);
  Fp(N+1) = fk(1);

  for j = 1:N/2
    tmp = circshift(Gp,-N/2+j-1);
    A(j, 1:N/2+j) = tmp(1:N/2+j);
  end
  A(N/2+1,1:end) = Gp(1:end);
  for j = 1:N/2
    A(N/2+1+j, j+1:end) = Gp(1:N+1-j);
  end

  Sk = A\Fp;
  sk = ifftshift(Sk(1:N))*N;
  sk = sk.';

end

