function out = petviashvili(Nin, nu, sigma, wavenumber, iter_max)
  global k N

  N = Nin;
  k = fftshift(-N/2:N/2-1);
  u = pi*(2*(0:N-1)/N - 1);
  
  % physical parameters and initial data
  a = 0.1;
  z = i*exp(-1.i*u).*(1 + a*exp(-1.i*wavenumber*u));
  %z = i*a*exp(-1.i*u);
  out = ifft(z);
  angular(z);
  
  dz = ifft(1.i*k.*fft(z)); 
  Nz = operator_N(z, nu, sigma);

  gamma = 1.0; 
  invK = 1./abs(k); invK(1) = 0.;
  %S = multiplier_S(z, nu, sigma)
  S2 = multiplier_S2(z, nu, sigma)
  for j = 1:iter_max
    z = (S2^gamma)*ifft(invK.*fft(operator_N(z, nu, sigma)));
    %zk = fft(z); zk(1) = 0; z = ifft(zk);
    %S = multiplier_S(z, nu, sigma)
    S2 = multiplier_S2(z, nu, sigma)
  end

  kw = wavenumber;
  out = z;
  figure(1)
  plot(u,real(z), u, imag(z));
  figure(2)
  plot(real(z), imag(z))
  figure(3)
  semilogy(k, abs(fft(z)/N))
end

function out = angular(in)
  global k 
  tmp = in.*conj(in);
  tmp = tmp.*ifft(abs(k).*fft(tmp));
  out = real(sum(tmp)/length(k));
end


function out = operator_N(in, nu, sigma) 
  % debugged and working properly
  global k N

  P   = 0.5*(1 - sign(k));
  dz  = fft(-1.i*k.*ifft(in));
  kz2 = ifft(abs(k).*fft(in.*conj(in))); 

  tmp1 = sigma*ifft( 1.i*k.*P.*fft(dz./abs(dz)));
  tmp2 = 0.5*nu^2*ifft(P.*fft(in.*kz2));

  %out = tmp1 + tmp2;
  
  filter_mask = P; filter_mask(N/2:N/2+N/8) = 0; 
  out = ifft(filter_mask.*fft(tmp1 + tmp2));
end

function out = multiplier_S(in, nu, sigma)
  global k N
 
  Nz = operator_N(in, nu, sigma);
  Mz = ifft(abs(k).*fft(in));

  num = real(sum(Mz(1:N).*conj(in(1:N))));
  den = real(sum(Nz(1:N).*conj(in(1:N))));
  
  out = abs(num/den); 

end


function out = multiplier_S2(in, nu, sigma)
  global k N

  P   = 0.5*(1 - sign(k));
  dz  = fft(-1.i*k.*ifft(in));
  kz  = ifft(abs(k).*fft(in));
  kz2 = ifft(abs(k).*fft(in.*conj(in))); 

  num = sum(in.*kz) + sigma*sum(abs(dz));
  den = 0.5.*nu^2*sum(abs(in).^2.*kz2);

  out = num/den;

end




