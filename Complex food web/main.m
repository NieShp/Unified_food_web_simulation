clear;
SIM = 1; RE=[];
for ff = 1:SIM
    % Plant and animal initial richness
    par.S_b = 20 ;   % initial plant richenss
    par.S_c = 40 ;   % initial animal richness
    
    % Holling type
    q = [1]; par.q  = q(randi([1, length(q)],1));
    
    % predator-prey mass ratio (PPMR)
    R_opt = 10.^2; par.R_opt  = R_opt(randi([1, length(R_opt)],1));
    
    % mineral (inorganic) nutrient supply concentration
    S0 = 10.^[-0.5:0.2:1.5] ; par.S  = S0(randi([1, length(S0)],1));
    
    % fraction of metabolic nutrient is egested organic detritus pool (indierct nutrient cycling)
    rho = 0.25; par.rho  = rho(randi([1, length(rho)],1));
    
    par.mu = 0.25;               % mineral nutrient supply rate - Brose 2008
    par.nu = 3*10^-4*365;        % mineral nutrient lose rate - Cherif & Loreau 2013
    par.nu_detritus = 8.4*10^-4*365;  % organic detritus nutrient lose rate - Cherif & Loreau 2013
    
    % species nutrient storage: max N/C ratio and min N/C ratio
    Q_max_animal = [0.2,0.3,0.4]; par.Q_max_species = Q_max_animal(randi([1, length(Q_max_animal)],1));   % maximal animal N:C ratio
    par.Q_max = [0.20*ones(par.S_b, 1); par.Q_max_species*ones(par.S_c, 1)] ;
    par.Q_min = [0.05*ones(par.S_b, 1); (par.Q_max_species-0.15)*ones(par.S_c, 1)] ;
    
    % N:C ratio of decomposer  (0.125 is from Manzoni & Porporato 2007)
    NC_ratio_decomposer = 0.125; par.NC_ratio_decomposer = NC_ratio_decomposer(randi([1, length(NC_ratio_decomposer)],1));
    
    warningState = warning('off', 'all');
     
    % predation and food structure
    par.sigma=2; par.ricker = 2; % width of the Ricker curve (higher value -&gt; narrower curve)
    par.range_b = [0 4];        % log_10 range; a larger range seems important to maintain more species (because it causes specialists of consumers)
    par.range_c = [2 10];       % log_10 range
    par.f_herbiv = 0;           % fraction of species that are strict herbivores
    par.f_pred = 0.00;          % fraction of species that are strict predators
    par.cutoff = 0.01;          % cutoff of the Ricker curve for setting a link between predator and prey
    % generate species body mass and food web topology based on body mass
    [par.mass, par.L] = foodweb_generate(par);
    % species group
    kP=find(sum(par.L,2)==0);    % plants
    kA=find(sum(par.L,2)~=0);    % animals
    kH=intersect(find(sum(par.L(:,kA),2)==0),kA); % Herbivores (only eat plants)
    kC=setdiff(kA,kH);                            % Carnivores (animal that can eat animals), including omnivores
    
    % importance between basal metabolism vs. stoichiometrically dominated metabolism for animal matabolism
    par.lambda = 0.5;
    
    % species metabolism rate
    par.m=zeros(par.S_b+par.S_c,1);
    par.m(1:par.S_b,1)= 0.138*par.mass(1:par.S_b).^(-0.25);  % Yoids & Innes 1992
    par.m(1+par.S_b:par.S_b+par.S_c,1)= 0.314*par.mass(1+par.S_b:par.S_b+par.S_c).^(-0.25); % Yodzis & Innes 1992
    
    % plant maximum growth rate
    par.r_max = 1 * par.mass(1:par.S_b).^(-0.25);   % Brown et al 2004
    
    % decomposer growth rate and dead rate from Cherif & Loreau 2013
    par.l = 0.3;
    par.x_M = 0.15;
    
    % plant nutrient taking up from mineral nutrient
    par.V =  zeros(par.S_b,1) + 1;
    par.K = rand(par.S_b,1) * 4 + 1;
    par.plant_affinity = par.V./par.K;    % plant nutrient affinity
    
    % feeding relationships (Rall et al. 2012; Schneider et al. 2016)
    % bij (attack rate), hij (handling time), c (predation interference) and assimilation
    % attack rate
    par.a = zeros(par.S_c+par.S_b);
    b0 = 0.45 ;
    beta_Cons=normrnd(0.47, 0.00,[par.S_b+par.S_c, 1]);
    beta_Prey=normrnd(0.15, 0.00,[par.S_b+par.S_c, 1]);
    par.a=b0*par.mass.^beta_Cons.*(par.mass.^beta_Prey)';
    % par.a(:,1:par.S_b)=b0*par.mass.^beta_Cons.*ones(1,par.S_b)*1;   % plant cannot move
    % handing time
    h0 = 0.0001;
    h_Cons=normrnd(-0.48, 0.00, [par.S_b+par.S_c, 1]);
    h_Prey=normrnd( 0.34, 0.00, [par.S_b+par.S_c, 1]);
    par.h=h0*par.mass.^h_Cons.*(par.mass.^h_Prey)';
    % predation interference
    par.c = 0.015 * ones(par.S_b+par.S_c, 1);

    
    
    % maximal assimilation rate
    par.e_max=zeros(par.S_b+par.S_c,par.S_b+par.S_c);
    par.e_max(:,1:par.S_b)= 0.85;       % when recource is plant
    par.e_max(:,1+par.S_b:end)= 0.85;   % wehn resource is animal
    
    % initial value of ODEs in food webs dynamics
    B0_carbon = 0 * rand(par.S_b+par.S_c,1) +  0.1;           % initial carbon stock
    B0_nutrient = B0_carbon.*(par.Q_max+par.Q_min)/2;         % initial nutrient stock
    BO_decomposer = 0.1*[1; par.NC_ratio_decomposer];         % initial decomposer community
    B0_detritus = 0.1*[1; par.NC_ratio_decomposer];           % initial detritus
    B0_mineral = 0.1;  % initial mineral nutrient
    B0 = [B0_carbon; B0_nutrient; BO_decomposer; B0_detritus; B0_mineral]; % initial value
    
    options = odeset( 'MaxStep', 1);
    [t_B,y_B]=ode15s(@(t,B) foodweb_dynmaics(B,par),[0:1:3*1e5], B0 );
    
    % calculate food web stcuture and ecosystem based on the last 30-steps
    Step = 30;
    clear primary_energy primary_nutrient Carbon_metabolism Nutrient_metabolism nutrient_excretion phi N_C_ratio Biomass_stock Nutrient_stock
    % Ecosystem properties calculation at equilibriums
    
    for step = 1:Step
        % calculation eccosystem functioning (e.g., feeding rates, productivity, realized assimialtion, etc.)
        [F, LL, primary_energy(step,:),primary_nutrient(step,:),Carbon_metabolism(step,:),...
            Nutrient_metabolism(step,:),nutrient_excretion(step,:),phi(step,:)] = Functioning_calculation(y_B(end-Step+step,:)',par);
        % biomass and nutrient stock
        Biomass = y_B(end-Step+step,1:par.S_b+par.S_c)';
        Nutrient = y_B(end-Step+step,1+par.S_b+par.S_c:2*(par.S_b+par.S_c))';
        N_C_ratio(step,:)=(Nutrient./Biomass)';
        Biomass(Biomass<1e-6)=0;Nutrient(Biomass<1e-6)=0;
        Biomass_stock(step,:) = [sum(Biomass(kP)),sum(Biomass(kH)),sum(Biomass(kC)),sum(Biomass)];
        Nutrient_stock(step,:) = [sum(Nutrient(kP)),sum(Nutrient(kH)),sum(Nutrient(kC)),sum(Nutrient)];
    end
    
    % detritus carbon and nutrient, as well as mineral nutrient
    decomposer_carbon = y_B(end-Step+1:end, 2*(par.S_b+par.S_c) + 1);     % mucrobial decomposer carbon stock (biomass)
    decomposer_nutrient = y_B(end-Step+1:end, 2*(par.S_b+par.S_c) + 2);   % mucrobial decomposer nutrient stock
    detritus_carbon = y_B(end-Step+1:end, 2*(par.S_b+par.S_c) + 3);       % detritus carbon stock
    detritus_nutrient = y_B(end-Step+1:end, 2*(par.S_b+par.S_c) + 4);     % detritus nutrient stock
    mineral_nutrient = y_B(end-Step+1:end, 2*(par.S_b+par.S_c) + 5);      % mineral nutrient
    
    % average ecosystem functioning
    Biomass_stock_mean = mean(Biomass_stock,1);
    Nutrient_stock_mean = mean(Nutrient_stock,1);
    primary_energy_mean = mean(primary_energy);
    primary_nutrient_mean = mean(primary_nutrient);
    Carbon_metabolism_mean = mean(Carbon_metabolism,1);
    Nutrient_metabolism_mean = mean(Nutrient_metabolism,1);
    nutrient_excretion_mean = mean(nutrient_excretion,1);
    decomposer_carbon_mean = mean(decomposer_carbon,1);
    decomposer_cnutrient_mean = mean(decomposer_nutrient,1);
    detritus_carbon_mean = mean(detritus_carbon,1);
    detritus_cnutrient_mean = mean(detritus_nutrient,1);
    mineral_nutrient_mesn = mean(mineral_nutrient,1);
    phi_mean = mean(phi,1);
    Functioning_mean = [Biomass_stock_mean,Nutrient_stock_mean,primary_energy_mean,primary_nutrient_mean,...
        Carbon_metabolism_mean,Nutrient_metabolism_mean,nutrient_excretion_mean,...
        decomposer_carbon_mean,decomposer_cnutrient_mean,...
        detritus_carbon_mean,detritus_cnutrient_mean,mineral_nutrient_mesn,phi_mean];
    
    % realized food web structure
    F=F.*Biomass; F(Biomass<1e-6,:)=[];F(:,Biomass<1e-6)=[];
    LL(Biomass<1e-6,:)=[];LL(:,Biomass<1e-6)=[];  % remove all extinct species
    % calculate tha survival plant and animal at the last step
    S_P=find(sum(LL,2)==0);    % survival plants
    S_A=find(sum(LL,2)~=0);    % survival animals
    S_H=intersect(find(sum(LL(:,S_A),2)==0),S_A);    % survival herbivores
    S_C=setdiff(S_A,S_H);                            % survival carnivores
    
    % average species realized N:C ratio
    N_C_ratio(:,Biomass<1e-6)=[];
    N_C_ratio_mean = mean(N_C_ratio,1);
    species_NC_ratio_mean = [mean(N_C_ratio_mean(S_P)),mean(N_C_ratio_mean(S_H)),mean(N_C_ratio_mean(S_C))];
    
    % survival species body mass
    Body_mass = par.mass; Body_mass(Biomass<1e-6)=[];
    plant_affinity = par.plant_affinity; plant_affinity(Biomass(1:par.S_b)<1e-6)=[];
    
   
    Struture = [length(S_P), length(S_H), length(S_C), size(LL,1)];
    
    % food web structure and functioning
    RE(ff,:) = [par.S, par.Q_max_species, par.rho, par.R_opt, par.q,...
        Struture, species_NC_ratio_mean, Functioning_mean];

end
clf
% plot for single simulation
if SIM == 1
    y_B(:,(1+par.S_b+par.S_c):2*(par.S_b+par.S_c)) = y_B(:,(1+par.S_b+par.S_c):2*(par.S_b+par.S_c))./y_B(:,1:(par.S_b+par.S_c));
    idx_p=find(y_B(end,(1:par.S_b))>1e-6);
    idx_a=find(y_B(end,(1+par.S_b):(par.S_b+par.S_c))>1e-6);    
    figure(1);
    subplot(231);plot(t_B, y_B(:,idx_p),'LineWidth',1.5);title('Plant carbon dynamics');set(gca,'FontName','Arial','FontSize',12,'LineWidth',1.5)
    subplot(232);plot(t_B, y_B(:,idx_p + par.S_b + par.S_c),'LineWidth',1.5);title('Plant N:C dynamics');set(gca,'FontName','Arial','FontSize',12,'LineWidth',1.5)
    subplot(234);plot(t_B, y_B(:,par.S_b + idx_a),'LineWidth',1.5);title('Amnimal carbon dynamics');set(gca,'FontName','Arial','FontSize',12,'LineWidth',1.5)
    subplot(235);plot(t_B, y_B(:,idx_a + par.S_b + par.S_b + par.S_c),'LineWidth',1.5);title('Animal N:C dynamics');set(gca,'FontName','Arial','FontSize',12,'LineWidth',1.5)
    subplot(233);plot(t_B, y_B(:,end-4:end-3),'LineWidth',1.5);title('Decomposer carbon & nutrient dynamics');set(gca,'FontName','Arial','FontSize',12,'LineWidth',1.5)
    subplot(236);plot(t_B, y_B(:,end-2:end),'LineWidth',1.5);title('Detritus mineral & nutrient');set(gca,'FontName','Arial','FontSize',12,'LineWidth',1.5)
    
end