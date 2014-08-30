clear 
clc;

global z h len dt A B C D E delta alp I II

% system parameters
m=1/2;
qr=0;
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
D=(Omega/2)*(i*epsilonp)/(deltac+i*kapa);
E=(Omega^2/4)/(deltac+i*kapa);
alp=sqrt(m*omegaz);
I=eye(len,len);
II=eye(2*len,2*len);

% differential parameters
Num=500;     % grid number
zstart=-10;
zend=10;
h=(zend-zstart)/Num;     % spacial step
z=zstart:h:zend-h;
len=length(z);
dt=-i/10^2;
profile on
[u,psi]=calculation;
profile off
profile viewer
