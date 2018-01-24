function [nlin_z,nlin_U] = DtUzeta_ETDRK4(zhat, U)

global etahat del2 mu KX KY Nx Ny F;

if abs(zhat(KX==0&KY==0))>1e-10
    error('zhat has power at (kx,ky)=(0,0)');
end

psihat = zhat./del2;

formstress = real(mean(conj(psihat(:)).*(1i*KX(:).*etahat(:)))) / (Nx*Ny);

q = +real(ifft2(zhat+etahat));
u = -real(ifft2(1i*KY.*psihat));
v = +real(ifft2(1i*KX.*psihat));

nlin_U = F - mu*U - formstress;
nlin_z = -1i*KX.*fft2((U + u).*q) - 1i*KY.*fft2(v.*q);

end