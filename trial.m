

psi_up=[1:0.1:50];
psi_down=[5.1:0.1:100];

P=D+E*h*conj(psi_down)*psi_up.';
Del=sqrt(delta^2+P*P');
W=exp(i*dt/2*Del);

%% 1st step evolution in position space
f=exp(-i*dt/2*C*z.^2).*(((delta-Del)^2*W+P*P'*W')*psi_up+P*(delta-Del)*(W-W')*psi_down);
g=exp(-i*dt/2*C*z.^2).*(P'*(delta-Del)*(W-W')*psi_up+((delta-Del)^2*W'+P*P'*W)*psi_down);

%% 2nd step evolution and start Fourier Transform
L=2*max(abs(z));
n=[-len/2:1:len/2-1];
k=2*n*pi/L;

% Fourier transform for f
cf=fftshift(fft(f));
cf=exp(-i*dt*(A*k.^2-B*k)).*cf;
F=ifft(fftshift(cf)); 

% Fourier transform for g
cg=fftshift(fft(g));
cg=exp(-i*dt*(A*k.^2+B*k)).*cg;
G=ifft(fftshift(cg));

%% 3rd step evolution and Find new psi_up and psi_down
npsi_up=(1/(2*Del*(Del-delta))^2)*exp(-i*dt/2*C*z.^2).*(((delta-Del)^2*W+P*P'*W')*F+P*(delta-Del)*(W-W')*G);
npsi_down=(1/(2*Del*(Del-delta))^2)*exp(-i*dt/2*C*z.^2).*(P'*(delta-Del)*(W-W')*F+((delta-Del)^2*W'+P*P'*W)*G);

psi=[npsi_up,npsi_down];