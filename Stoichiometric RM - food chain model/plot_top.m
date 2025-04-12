clear
clf

par.a = 0.81;  % attack rate
par.h = 0.1;  % handlign time
par.r= 1.2;    % resource growth rate
par.m= 0.27;   % consumer mortality
par.Q0 = 0.0038;  % resource minimal N:C ratio
par.e=0.8;     % consumer maximal assimialtion rate

a = par.a; h= par.h; x_H = par.m; Q_Pmin = par.Q0; e = par.e; r = par.r;
color1 = [144/255, 238/255, 144/255];  % #70AD47
color2 = [91/255 156/255 213/255];     % #5B9BD5

for u = 1 : 4
    
    if u ==1; N_tot = 10^-1.6; Q_H = 0.045; subplot(2,4,1);  end
    if u ==2; N_tot = 10^-1.1; Q_H = 0.045; subplot(2,4,2);   end
    if u ==3; N_tot = 10^-1.6; Q_H = 0.01;  subplot(2,4,3);  end
    if u ==4; N_tot = 10^-1.1; Q_H = 0.01;  subplot(2,4,4);   end
    
    
    A = -a*h*Q_Pmin / N_tot;
    B = -a*h*Q_H / N_tot;
    C = a*Q_H / (r*N_tot);
    D = (a*h*N_tot - Q_Pmin) / N_tot;
    E = -(r*Q_H + a*N_tot) / (r*N_tot);
    

    x = x_H / (a*(e - h*x_H));
    % herbivore zero growth isocline
    C_P2 = linspace(x, e*N_tot/(h*x_H * Q_H) - 1/(a*h) , 1e5);
    C_H2 = -h*x_H/e * C_P2 + N_tot / Q_H - x_H/(e*a);
    

    x_intercept = N_tot/Q_H;  
    y_intercept = N_tot/Q_H;  
    fill([0, 0, x_intercept], [N_tot / Q_H, 0, 0],[0.9, 0.9, 0.9]); hold on; % 淡蓝色填充
    C_P2(C_H2<0)=[];C_H2(C_H2<0)=[];
    plot(C_P2, C_H2,  'Color', color2,  'LineWidth', 2); hold on
    yLimits = ylim; line([x x], [0 -h*x_H/e * x + N_tot / Q_H - x_H/(e*a)], 'Color', color2, 'LineWidth', 2); hold on
    
    
    %plant zero growth isocline
    C_P = linspace(0, N_tot/Q_Pmin , 1e5);
    A1 = C;
    B1 = B * C_P + E;
    C1 = A * C_P.^2 + D * C_P + 1;
    C_H1 = (- B1 - sqrt(B1.^2 - 4*A1 .* C1)) ./ (2*A1);
    plot(C_P, C_H1, 'Color', color1, 'LineWidth', 2); hold on;
    
    

    C_P2 = linspace(0, N_tot/Q_H , 1e3);
    C_H2 = N_tot / Q_H - C_P2;
    plot(C_P2, C_H2,  ':', 'Color', [1, 0.6471, 0], 'LineWidth', 2); hold on
    
    
    
    set(gca, 'FontName', 'Arial', 'FontSize', 14);

    
    if u == 1; xlim([0 6]);  end
    if u == 2; xlim([0 4]); end
    if u == 3; xlim([0 4]);  end
    if u == 4; xlim([0 4]);   end
    
end

