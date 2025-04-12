function [mass, LL]= foodweb_generate(par)
% generate foodweb structure based on the "body-mass algorithm"

% "mass" - a positive vector with (S_b + S_c) elements
%       it includes the body mass of basal species (first S_b elements) and
%       consumers (remaining S_c elements);
%       Note: the first S_b values (for basal) are in increasing order;
%       same for the remaining S_c value
% "LL" - a nonnegative matrix of (S_b + S_c)*(S_b + S_c)
%       where L(i,j) is the feeding probability of species i on j.
%       Thus, the first S_b rows (i.e. basal species) are all 0
%       Values smaller than cutoff are set to 0, which results in a
%       lower-triangle matrix (recalling that body mass have been sorted).

flag = 0;   % in case the generated foodweb contains uncontrolled basal species and consumers with no prey

while flag < 1
    
    % generate biomass for basal (mass_b) and consumer species (mass_c)
    % Note: in case some values are out of interest, first generate 0.5 times more species,
    % and then select S_b and S_c, respectively)
    
    tmp_b = par.range_b(1) + rand(par.S_b,1).*(par.range_b(2)-par.range_b(1));
    mass_b = sort(10.^tmp_b);
    tmp_c = par.range_c(1) + rand(par.S_c,1).*(par.range_c(2)-par.range_c(1));
    mass_c = sort(10.^tmp_c);
    
    % mass_b(1) = 10^par.range_b(1); mass_c(1) = 10^par.range_c(1);
    mass = [mass_b; mass_c];
    % calculate the feeding probability based on Ricker function
    LL = zeros(par.S_b+par.S_c);
    for k = 1:par.S_c
        LL(k+par.S_b,:) = (mass_c(k)/par.R_opt ./mass .* exp(1-mass_c(k)/par.R_opt./mass)).^par.ricker;
        
        tmp = rand(1);
        if(tmp < par.f_herbiv)
            LL(k+par.S_b,(par.S_b+1):(par.S_b+par.S_c)) = 0;
        elseif(tmp < par.f_herbiv+par.f_pred)
            LL(k+par.S_b,1:par.S_b) = 0;
        end
    end
    LL(LL<par.cutoff) = 0;
    flag = all(sum(LL((par.S_b+1):(par.S_b+par.S_c),:),2) > 0) & all (sum(LL(:,1:par.S_b),1) > 0);
end
end

