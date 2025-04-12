clear
par.a = 0.81;     % attack rate
par.h = 0.1;      % handlign time
par.r= 1.2;       % resource growth rate
par.m= 0.27;      % consumer mortality
par.Q0 = 0.0038;  % resource minimal N:C ratio
par.e=0.8;        % consumer maximal assimialtion rate

R = [];

for i = 1:4
    if i == 1;  NT = -1.6; q = 0.045;  B0 = [0.002; 0.1]; subplot(2,4,5); end  % (a)
    if i == 2;  NT = -1.1; q = 0.045;  B0 = [0.8; 1.46]; subplot(2,4,6);  end  % (b)
    if i == 3;  NT = -1.6; q = 0.01;   B0 = [0.8; 0.5];  subplot(2,4,7);  end  % (c)
    if i == 4;  NT = -1.1; q = 0.01;   B0 = [0.1; 5.88]; subplot(2,4,8);  end  % (d)
    
    
    par.q = q;           % consumer N:C ratio
    par.N = 10.^NT;      % total nutrient

    
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'MaxStep', 0.1);    
    if i ==1; Time = 300; [t_B, y_B] = ode45(@(t, B) SimpleFoodChain(B, par), [0, Time], B0,options); end
    if i ==2; Time = 300; [t_B, y_B] = ode45(@(t, B) SimpleFoodChain(B, par), [0, Time], B0,options); end
    if i ==3; Time = 300; [t_B, y_B] = ode45(@(t, B) SimpleFoodChain(B, par), [0, Time], B0,options); end
    if i ==4; Time = 200; [t_B, y_B] = ode45(@(t, B) SimpleFoodChain(B, par), [0, Time], B0,options); end
    RE = [t_B, y_B];
    
    if i == 1; RE(:,4)=1; end  % (a)
    if i == 2; RE(:,4)=2; end  % (b)
    if i == 3; RE(:,4)=3; end  % (c)
    if i == 4; RE(:,4)=4; end  % (d)
    
    R = [R; RE];
    clear RE
    

    color1 = [144/255, 238/255, 144/255];  % #70AD47
    color2 = [91/255 156/255 213/255];     % #5B9BD5
    plot(t_B, y_B(:,1), 'Color', color1, 'LineWidth', 2); hold on
    plot(t_B, y_B(:,2), 'Color', color2, 'LineWidth', 2); hold on
    if i == 2
        legend({'Plant', 'Herbivore'})
    end
    if i == 4
        ylim([0 8])
    end
    set(gca, 'FontName', 'Arial', 'FontSize', 14);

end


set(gcf, 'Position', [100, 100, 1400, 600]); % 图形尺寸为1200x900


