function H=hamiltonian(psi_up,psi_down)

global z h len dt A B C D E delta I

H1=diag(2*A/h^2+delta+C*z.^2)+(-A/h^2+i*B/(2*h))*diag(ones(len-1,1),-1)+(-A/h^2-i*B/(2*h))*diag(ones(len-1,1),1);
H2=diag(2*A/h^2-delta+C*z.^2)+(-A/h^2-i*B/(2*h))*diag(ones(len-1,1),-1)+(-A/h^2+i*B/(2*h))*diag(ones(len-1,1),1);
P1=(D+E*h*conj(psi_down)*psi_up.')*I;
P2=conj(P1);

H=[H1,P1;P2,H2];

end
