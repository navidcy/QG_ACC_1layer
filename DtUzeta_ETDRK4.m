function [nlin_z,nlin_U] = DtUzeta_ETDRK4(zhat, U)

global etahat detadx detady del2 mu KX KY Nx Ny F;

if abs(zhat(KX==0&KY==0))>1e-10
    error('zhat has power at (kx,ky)=(0,0)');
end

psihat = zhat./del2;

formstress = real(mean(conj(psihat(:)).*(1i*KX(:).*etahat(:)))) / (Nx*Ny);

dpsidx  = real(ifft2(1i*KX.*psihat));
dpsidy  = real(ifft2(1i*KY.*psihat));
dzetadx = real(ifft2(1i*KX.*zhat));
dzetady = real(ifft2(1i*KY.*zhat));

J  = dpsidx.*(dzetady+detady) - dpsidy.*(dzetadx+detadx);

nlin_U = F - mu*U - formstress;
nlin_z = fft2( -J -U*(dzetadx+detadx) );

end