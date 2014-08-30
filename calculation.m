function [u,psi]=calculation(alp,z,h)
len = length(z);
%% Exact ground state
exact_psi=1/sqrt(2)*sqrt(alp)*exp(-alp^2*(z).^2/2)/pi^(1/4);

%% Initial psi_up and psi_down

% initial psi_up
psi_up=sqrt(alp)*exp(-alp^2*(z-1).^2/2)/pi^(1/4);
% initial psi_down
psi_down=sqrt(alp)*exp(-alp^2*(z+1).^2/2)/pi^(1/4);

% normalization
N=sqrt(h*sum(abs(psi_up).^2)+h*sum(abs(psi_down).^2));
psi_up=psi_up/N;
psi_down=psi_down/N;
psi=[psi_up,psi_down];

% initial chemical potential
mu=h*conj([psi_up,psi_down])*hamiltonian(psi_up,psi_down)*[psi_up,psi_down].';

%% Find stable solution

d=1;
ncount=1;

while d>10^(-10)
    
    % new psi after evolution
    psi=evolution(psi_up,psi_down);
    
    % Normalization
    N=sqrt(h*sum(abs(psi).^2));
    psi=psi/N;
    
    % new psi_up and psi_down after evolution
    psi_up=psi(1:len);
    psi_down=psi(len+1:2*len);
    
    % chemical potential
    mup=h*conj(psi)*hamiltonian(psi_up,psi_down)*psi.';
    d=abs(mu-mup);
    mu=mup;
    
    % draw psi_up and psi_down
    if mod(ncount,5)==0
        plot(z,abs(psi_up).^2,'r',z,abs(psi_down).^2,'b',z,abs(exact_psi).^2,'y')
        drawnow
    end
     ncount=ncount+1;
       
end

u=mu

end