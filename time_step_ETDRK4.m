function  [zhatnew,Unew] = time_step_ETDRK4(zhat,U,eL,eL2,Q,fu,fab,fc,dt);

[nlin0_z,nlin0_U] = DtUzeta_EDTRK4(zhat,U);
k1z = eL2.*zhat + Q.*nlin0_z;   k1U = nlin0_U;

[nlin1_z,nlin1_U] = DtUzeta_EDTRK4(k1z,U+k1U*dt/2);
k2z = eL2.*zhat + Q.*nlin1_z;   k2U = nlin1_U;

[nlin2_z,nlin2_U] = DtUzeta_EDTRK4(k2z,U+k2U*dt/2);
k3z = eL2.*k1z + Q.*(2*nlin2_z-nlin0_z);   k3U = nlin2_U;

[nlin3_z,nlin3_U] = DtUzeta_EDTRK4(k3z,U+k3U*dt);
                                   k4U = nlin3_U;

zhatnew  = eL.*zhat + fu.*nlin0_z + 2*fab.*(nlin1_z+nlin2_z) +fc.*nlin3_z ;
Unew     = U + (k1U+2*k2U+2*k3U+k4U)*dt/6;

end