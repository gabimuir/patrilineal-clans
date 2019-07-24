% SCUDEM Patrilineal Clans Code
% Oct 2018
% Gabi, Nguyen, Heather

%Initial Conditions
%Tribe 1
M1 = 50;
F1 = 50;
%Tribe 2
M2 = 25;
F2 = 25;

% Rate Constants
birth_rate = 0.01; %from 0.01
dampening_rate = 0.005;
max = 300;
marriage_rate_1 = 0.001; %c1
marriage_rate_2 = 0.001; %c2
combat_rate_1 = 0.003; %q
combat_rate_2 = 0.003;

%time scale & options
MaxTime = 50; % time in [units]
tSpan = [0:0.1:MaxTime]; % [start: stepsize: end]
Options1 = odeset('NonNegative', 1, 'RelTol',1e-4,'AbsTol',1e-5);

figure
k = 4;
j = linspace(0.001, 0.01, k);
for subp = 1:k
    combat_rate_1 = j(subp);
    
    InitVals = [M1, F1, M2, F2];
    constants = [birth_rate, dampening_rate, max, marriage_rate_1, marriage_rate_2, combat_rate_1, combat_rate_2];

    
    %run
    [t, Solutions] = ode15s( @(t,y)popdym(t,y,constants), tSpan, InitVals, Options1);

    %plot
    subplot (1,k,subp);
    plot(tSpan, Solutions(:,1),'LineWidth', 2);  %M1
    hold on
    plot(tSpan, Solutions(:,2),'LineWidth', 2);  %F1
    hold on
    plot(tSpan, Solutions(:,3),'LineWidth', 2);  %M2
    hold on
    plot(tSpan, Solutions(:,4), 'LineWidth', 2);  %F2
    hold on
    ylim([0 max/2.5])
    xlim([0 MaxTime])
    xlabel('Time')
    ylabel('Population')
    titlestr = sprintf('q1 = %s', num2str(combat_rate_1));
    title([ titlestr ])
    if subp == 1
        legend('Males Tribe 1', 'Females Tribe 1', 'Males Tribe 2', 'Females Tribe 2')
    end
end

suptitle('Varying combat morbidity for clan 1. (q2 = 0.003)')
% string = sprintf('IC: M1=%s F1=%s M2=%s F2=%s Parameters: r1=%s, r2=%s, max=%s, c1=%s, c2=%s, q=%s', num2str(M1), num2str(F1),num2str(M2),num2str(F2),num2str(birth_rate), num2str(dampening_rate), num2str(max), num2str(marriage_rate_1),num2str(marriage_rate_2), "varying");
% suptitle(string)

%Returns the change in population each step
function [dPdt] = popdym(t, current_vals, consts)
    %Get the current_vals
    c = num2cell(current_vals); %, ones(size(current_vals,1),1), ones(1,size(current_vals,2)));
    [M1, F1, M2, F2] = c{:}; %the same order as the initial conditions

    birth_rate = consts(1);
    dampening_rate = consts(2);
    max = consts(3);
    marriage_rate_1 = consts(4); 
    marriage_rate_2 = consts(5);
    combat_rate_1 = consts(6);
    combat_rate_2 = consts(7);
        
    %Birth and death population dynamics
    dF1b = birth_rate*F1*M1*(1 - F1/(max-M1)) - dampening_rate*F1*F1; 
    dF2b = birth_rate*F2*M2*(1 - F2/(max-M2)) - dampening_rate*F2*F2;
    dM1b = birth_rate*F1*M1*(1 - M1/(max-F1)) - dampening_rate*M1*M1;
    dM2b = birth_rate*F2*M2*(1 - M2/(max-F2)) - dampening_rate*M2*M2;
    
    %combat
    dcombat1 = combat_rate_1*M1*M2;
    dcombat2 = combat_rate_2*M1*M2;
    
    %marriage (outside)
    dF1M2m = marriage_rate_1 * F1*M2; % based just on the female interactions with their males
    dF2M1m = marriage_rate_2 * F2*M1;
    
    dPdt = [
        dM1b - dcombat1; % M1
        dF1b - dF1M2m + dF2M1m; % F1
        dM2b - dcombat2; % M2
        dF2b + dF1M2m - dF2M1m %F2
        ];
end

