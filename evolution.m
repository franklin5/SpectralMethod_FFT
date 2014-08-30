function psi=evolution(psi_up,psi_down)

global z len dt A B C  delta

%% 1st step evolution in position space
f=exp(-1i*dt/2*(C*z.^2+delta)).*psi_up;
g=exp(-1i*dt/2*(C*z.^2-delta)).*psi_down;

%% 2nd step evolution and start Fourier Transform
L=2*max(abs(z));
n=-len/2:1:len/2-1;
k=2*n*pi/L;

% Fourier transform for f
cf=fftshift(fft(f));
cf=exp(-1i*dt*(A*k.^2+B*k)).*cf;
F=ifft(fftshift(cf)); 

% Fourier transform for g
cg=fftshift(fft(g));
cg=exp(-1i*dt*(A*k.^2-B*k)).*cg;
G=ifft(fftshift(cg));

%% 3rd step evolution and Find new psi_up and psi_down
npsi_up=exp(-1i*dt/2*(C*z.^2+delta)).*F;
npsi_down=exp(-1i*dt/2*(C*z.^2-delta)).*G;

psi=[npsi_up,npsi_down];

end
