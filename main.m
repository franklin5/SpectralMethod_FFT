clear
clc
close all
global z h len dt A B C D E delta alp I1 II

% system parameters
m=1; % why it was 1/2 in the past?
qr=1;
omegaz=1;
delta=0;
Omega=0;
epsilonp=1;
deltac=0;
kapa=1;

% parameters in calculation
A=1/(2*m);
B=qr/m;
C=(1/2)*m*omegaz^2;
D=(Omega/2)*(1i*epsilonp)/(deltac+1i*kapa);
E=(Omega^2/4)/(deltac+1i*kapa);
alp=sqrt(m*omegaz);

% differential parameters
Num=2^8;     % grid number
zstart=-10;
zend=-zstart;
h=(zend-zstart)/Num;     % spacial step
z=zstart:h:zend-h;
len=length(z);
dt=-1i/10^2;
I1=eye(len,len);
II=eye(2*len,2*len);
%%
%profile on
[u,psi]=calculation(alp,z,h);
%profile off
%profile viewer
