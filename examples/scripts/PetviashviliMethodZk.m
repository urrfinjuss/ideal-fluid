function [out] = PetviashviliMethodZk(Nin, Omega_in, A, k0, nsteps)
% this code doesn't work, and I suspect this is 
% because y_u/|z_u| is computed not very accurately.
% Two things:  -- do iterations in k-space
%	       -- write z_u/|z_u| with Fourier series
% N = 16384; z = PetviashviliMethod(N, 1.06, 0.05, 4, 1000);
% N = 65536; z = PetviashviliMethod(N, 1.2, 0.25, 4, 4000);  
global k Ak Pk Hk Omega Sigma iLk



  graphics_toolkit("gnuplot")
  k = fftshift(-Nin/2:Nin/2-1);
  u = pi*(2*(0:Nin-1)/Nin - 1);
  Ak = -1i./k;
  Ak(1) = 0;
 
  Pk = 0.5*(1-sign(k));
  Hk = 1.i*sign(k);
  iLk = 1./(1-abs(k)); iLk(2) = 0; iLk(end) = 0;
  
  Sigma = 1;
  Omega = Omega_in;
  zk = zeros(1,Nin);
  zk(end) = -1i;
  zk(end-k0) = 1i*A;
  
  for j = 1:nsteps

    S = (zk*zk')/(nonlinear(zk)*zk');   
    zk = S*nonlinear(zk);

    %mu  = 2.*pi*sum(abs(k).*abs(yk).^2);
    %yk = S.*yk;
    %zk = -1.*(Hk - 1.i).*yk;
    %yk = yk/sqrt(mu);
    %abZ2 = abs(ifft(zk)*Nin).^2; 
    %mu  = 2.*pi*sum(abs(k).*abs(yk).^2);
    %per = 2*pi*real(sum(fft(abs(ifft(zk)))));
    %mJ  = 2.*pi*sum(abs(k).*abs(fft(abZ2/Nin)).^2);
    % self-consistent lambda
    %Lambda = 0.5*(per - 0.5*Omega^2*mJ)/mu;
    
    %zk = -1.*(Hk - 1.i).*yk;
    %mu = 0.5*sum(abs(k).*abs(zk).^2);
    %fprintf('Petviashvili iteration %4d: S = %.12e\tmu = %.12e\tPerimeter = %.12e\tmJ = %.12e\n', j, S, mu, per, mJ);
    fprintf('Petviashvili iteration %4d: S = %.12e\n', j, S);
  end
  res = zk - nonlinear(zk);
  


  NormR = sqrt(sum(abs(res).^2));
  printf('Norm of Residual = %.12e\n', NormR);
  
  out = ifft(zk)*Nin;
  % prepare wave for the dynamic code:
  dz = ifft(1i*k.*fft(out));
  dPhi = 1i*Omega.^2*ifft(1i*k.*Pk.*fft(abs(out).^2));
  
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


function outk = nonlinear(ink)
 global k Ak Pk Hk Omega

 N = length(k);
 z = ifft(ink)*N;
 dz = ifft(1i*k.*ink)*N;
 kz2 = ifft(abs(k).*fft(abs(z).^2));

 %outk = 1i*Pk.*fft(exp(1i.*arg(dz)))/N;
 Tmp = ratio_fourier( ifftshift(fft(dz)), ifftshift(fft(abs(dz))));
 outk = 1i*Pk.*Tmp/N;
 outk = outk + 0.5i*Omega^2*Pk.*Ak.*fft(z.*kz2)/N;
 outk(N/4:3*N/4) = 0; 

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

