function psi=evolution(psi_up,psi_down)

global z h len dt A B C D E delta

%% 1st step evolution in position space
f=exp(-i*dt/2*(C*z.^2+delta)).*psi_up;
g=exp(-i*dt/2*(C*z.^2-delta)).*psi_down;

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
npsi_up=exp(-i*dt/2*(C*z.^2+delta)).*F;
npsi_down=exp(-i*dt/2*(C*z.^2-delta)).*G;

psi=[npsi_up,npsi_down];


%{
L=2*max(abs(x));
N=Num;
n = [-N/2:1:N/2-1].';           % Indices
e = -4*n.*n*pi*pi/L/L;         % Squares of wavenumbers.
c0=c;
u=psi.';
V=V.';

u = exp(-dt*1i*(V+c0*(abs(u).*abs(u)))/2.0).*u; % propagate in physical space
cu = fftshift(fft(u));                 % Take Fourier transform
cu = exp(dt*1i*e/2).*cu;                   % Advance in Fourier space
u = ifft(fftshift(cu));                % Return to physical space 
u = exp(-dt*1i*(V+c0*(abs(u).*abs(u)))/2.0).*u; % propagate in physical space                                        

% u=1/sqrt(sum(abs(u).^2*h))*u;
u=u.';
%}

end
