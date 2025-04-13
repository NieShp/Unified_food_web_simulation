function [F, LL, primary_energy,primary_nutrient,Carbon_metabolism,Nutrient_metabolism,nutrient_excretion,phi] = Functioning_calculation(y,par)
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

% Type I
% phi = par.l * (N_D - N_M/C_M * C_D);

% Type II
phi = par.l * (N_D * C_M - C_D * N_M);

% metabolism
kP=find(sum(LL,2)==0);    % plants
kA=find(sum(LL,2)~=0);    % animals
kH=intersect(find(sum(LL(:,kA),2)==0),kA);    % Herbivores (only eat plants)
kC=setdiff(kA,kH);                            % Carnivores (animal that can eat animals), including omnivores
mc = [mc_P; mc_A]; mn = [mn_P; mn_A]; 
Carbon_metabolism =[sum(C(kP).*mc(kP)),sum(C(kH).*mc(kH)),sum(C(kC).*mc(kC)),sum(mc.*C)];
Nutrient_metabolism =[sum(Q(kP).*C(kP).*mn(kP)),sum(Q(kH).*C(kH).*mn(kH)),sum(Q(kC).*C(kC).*mn(kC)),sum(Q.*C.*mn)];

% Gross primary productivity 
primary_energy = sum(par.r_max.*G.*C(indx_P));
primary_nutrient = sum(s.*C(indx_P));

% nutrient excretion (unassimilated nutrient)
eff_N(1:n1,:)=[]; 
nutrient_excretion = sum(sum((1-eff_N).*Q'.*F(indx_A,:).*C(indx_A), 2));





