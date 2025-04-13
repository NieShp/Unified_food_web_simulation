function dydt=foodweb_dynmaics(y,par)
% from top to the bottom
% Carbon stock of plant, animal; nutrient stock of plant, animal; 
% Carbon and nutrient of stock microbial decomposer; 
% Carbon and nutrient of stock detritus; 
% Nutrient stock of mineral nutrient

n1=par.S_b;    % plant diversity
n2=par.S_c;    % animal diversity
indx_P = 1:n1; % plant species index 
indx_A = (1+n1):n1+n2;    % animal species index 
C=y(1:n1+n2);             % species carbon stock (biomass)
N=y(1+n1+n2:2*(n1+n2));   % species nutrient carbon 
C_M = y(2*(n1+n2) + 1);   % mucrobial decomposer carbon stock (biomass)
N_M = y(2*(n1+n2) + 2);   % mucrobial decomposer nutrient stock
C_D = y(2*(n1+n2) + 3);   % detritus carbon stock
N_D = y(2*(n1+n2) + 4);   % detritus nutrient stock
L = y(2*(n1+n2) + 5);     % mineral nutrient stock
Q=N./C;                   % Species Nutrient: Carbon ratio
Q(intersect(find(Q>par.Q_max),find(C>1e-6))) = par.Q_max(intersect(find(Q>par.Q_max),find(C>1e-6)));
Q(intersect(find(Q<par.Q_min),find(C>1e-6))) = par.Q_min(intersect(find(Q<par.Q_min),find(C>1e-6)));
C(C<1e-6)=0;N = Q.*C;  % if species biomass is less than 1e-6, is considered as extinct
C_M(C_M<1e-6)=0;       % if decomposer biomass is less than 1e-6, is considered as extinct

% feeding rate (F: functioning response)
LL = par.L;             % food webs topological structure
a=par.a; a(LL==0)=0;    % attack rate
h=par.h; h(LL==0)=0;    % handling time
c=par.c;                % predation interference
F = (a.*real(C'.^par.q))./(1+c.*C+sum(a.*h.*real(C'.^par.q),2))./par.mass;  % feeding rates

% species metabolism rate (carbon & nutrient)
mc_P=par.m(indx_P).*((Q(indx_P)./par.Q_min(indx_P)));   % plant C metabolism
mn_P=mc_P;                                              % plant N metabolism
mc_A=par.m(indx_A).*(((Q(indx_A)-par.Q_min(indx_A))+par.lambda * (par.Q_max(indx_A) - Q(indx_A)))./(par.Q_max(indx_A)-par.Q_min(indx_A))); % animal C metabolism
mn_A=par.m(indx_A).*(((Q(indx_A)-par.Q_min(indx_A))+par.lambda * (Q(indx_A) - par.Q_min(indx_A)))./(par.Q_max(indx_A)-par.Q_min(indx_A))); % animal N metabolism

% plant nutrient dependent carbon growth factor
G =(1-par.Q_min(indx_P)./Q(indx_P));
% plant nutrient taking up rate
s =(par.V.*L./(par.K+L)).*(1-Q(indx_P)./par.Q_max(indx_P));

% food assimilation rate for C & N
eff_C = par.e_max .* min(1, (Q'./Q).* ((Q-par.Q_min)./(par.Q_max-par.Q_min)));
eff_N = par.e_max .* min(1, (Q'./Q).* ((par.Q_max-Q)./(par.Q_max-par.Q_min))); 

eff_C(1:n1,:)=[];eff_N(1:n1,:)=[];

% overall food webs danamics
dCPdt = par.r_max.*G.*C(indx_P) - mc_P.*C(indx_P) - sum(F(:, indx_P).*C, 1)' ;
dCAdt = sum(eff_C.*F(indx_A, :).*C(indx_A), 2) - mc_A.*C(indx_A) - sum(F(:,indx_A).*C, 1)' ;

dNPdt = s.*C(indx_P) - sum(F(:, indx_P).*Q(indx_P)'.*C, 1)' - mn_P.*N(indx_P) ;
dNAdt = sum(eff_N.*Q'.*F(indx_A,:).*C(indx_A), 2) - sum(F(:,indx_A).*Q(indx_A)'.*C, 1)' - N(indx_A).*mn_A ;

% Type I
% phi = par.l * (N_D - N_M/C_M * C_D);
% dC_Mdt = par.l * C_D - par.x_M * C_M;
% dN_Mdt = par.l * N_D - phi - par.x_M * N_M;

% Type II
phi = par.l * (N_D * C_M - C_D * N_M);
dC_Mdt = par.l * C_D * C_M - par.x_M * C_M ;
dN_Mdt = par.l * N_D * C_M - phi - par.x_M * N_M ;


dC_Ddt = sum(sum((1-eff_C).*F(indx_A,:).*C(indx_A), 2))  -  par.l * C_D * C_M - par.nu_detritus*C_D;
dN_Ddt = sum(sum((1-eff_N).*Q'.*F(indx_A,:).*C(indx_A), 2)) + par.rho*(sum(mn_P.*N(indx_P))+sum(N(indx_A).*mn_A))...
    - par.l * N_D * C_M - par.nu_detritus*N_D;

A = par.mu*(par.S - L);
if A <0
    A = 0;
end
dLdt = A + phi + (1-par.rho)*(sum(mn_P.*N(indx_P))+sum(N(indx_A).*mn_A)) - sum(sum(s.*C(indx_P))) - par.nu*L;

% species whose biomass (carbon stock) is less than 1e-6 is considered as extinct
dCPdt(C(indx_P)<1e-6)=0;dNPdt(C(indx_P)<1e-6)=0;
dCAdt(C(indx_A)<1e-6)=0;dNAdt(C(indx_A)<1e-6)=0;
dC_Mdt(C_M<1e-6)=0;dN_Mdt(C_M<1e-6)=0;
dydt = real([dCPdt;dCAdt;dNPdt;dNAdt;dC_Mdt;dN_Mdt;dC_Ddt;dN_Ddt;dLdt]);

end