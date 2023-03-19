%%
close all; clear; clc

%% Procedure
% 1st step: Define numbers of stations and slots in the simulation
N_stations = 20;
N_slots = 10000; 
% 2nd step: Predefine packet arrival (lambda)
packet_arrival = [1:20 25:5:50 60:10:200];
% 3rd step: Define variables to save results
G = zeros(1,length(packet_arrival));
S_analyt = zeros(1,length(packet_arrival));
S_sim = zeros(1,length(packet_arrival));

% 4th step: Select current lambda value used for Poisson
for lambda = 1:length(packet_arrival)

    for slot = 1:N_slots
        % 5th step: Generate initial packet arrival for each station (only for first slot)
        if slot == 1
            StationTimeToTx = poissrnd(1,1,N_stations);
        end
        % 6th step: Determine if there are stations that are supposed to transmit
        g = G(lambda);
        for station = 1:N_stations
            if StationTimeToTx(station) == 0
                G(lambda) = G(lambda) + 1;
            end
        end
        g = G(lambda) - g;
        if g == 1
            S_sim(lambda) = S_sim(lambda) + 1;
        end
        % 7th step: Check for collisions
        for station = 1:N_stations
            if StationTimeToTx(station) == 0
                StationTimeToTx(station) = poissrnd(packet_arrival(lambda));
            else
                StationTimeToTx(station) = StationTimeToTx(station) - 1;
            end
        end
    end

    % 8th step: Save results for current lambda %
    G(lambda) = G(lambda)/N_slots;
    S_sim(lambda) = S_sim(lambda)/N_slots;
    S_analyt(lambda) = G(lambda)*exp(-G(lambda));
end

probability = 100*S_sim./G;
probability_complement = 100*ones(1,length(probability)) - probability;

% 9th step: Plot the results
LS=12;
figure(1);
plot(G,S_analyt,'b');
hold on
plot(G,S_sim,'r');
legend('Slotted ALOHA: analyt','Slotted ALOHA: simulated','Location','NorthEast');
xlabel('G [-]');
ylabel('S [-]');
xlim([0 10])
ylim([0 0.5])
set(gca,'FontWeight','normal',...
    'FontSize',LS);%12
grid on;
print('-f1', '-r300', '-dpng','Results')

LS=12;
figure(2);
plot(G,probability_complement,'b');
hold on
plot(G,probability,'r');
legend('Collision','Succesful transmission','Location','NorthEast');
xlabel('G [-]');
ylabel('Probability [%]');
xlim([0 10])
ylim([0 100])
set(gca,'FontWeight','normal',...
    'FontSize',LS);%12
grid on;
print('-f1', '-r300', '-dpng','Results')
