%%
close all; clear; clc

%% Flooding routing algorithm with improvement and gossiping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. Step: Define basic parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Simulation_area = [500 500];  % Simulation area in m 
No_nodes = 500; % Definition of number of nodes
Sim_drops = 1000;  % Definition of number of drops (for outer loop)
Max_comm_d = 10:5:100; % Definition of max. comm. distance in m (for inner loop) 

p = 0.80;    % Gossiping factor p (probability of message forwarding)
TTL = 100;  % Maximum number of hops before the message is discarded 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. Step: Define variables to save results%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Total_messages = zeros(19,3);
Total_hops = zeros(19,3);

% START OF OUTER LOOP (For averaging over simmulation drops)
%---------------------------------------------------
for simulation_number = 1:Sim_drops
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %3. Step: Randomly generate positions of each node%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nodes_position = zeros(No_nodes,2);
    for node=1:No_nodes
        Nodes_position(node,1) = rand*Simulation_area(1,1); % x coordinates
        Nodes_position(node,2) = rand*Simulation_area(1,2); % y coordinates
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %4. Step: Select randomly one source (S) and one destination (D)%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Index_S_D = zeros(2,1); % parameter to store indexes of S and D
    Index_S_D(1,1) = ceil(rand*No_nodes); % Index of S node
    D = 0;
    while D == 0
        Index_S_D(2,1) = ceil(rand*No_nodes); % Index of D node
        if Index_S_D(2,1) ~= Index_S_D(1,1) % Check if S and D are different   
            D = 1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %5. Step: Calculation of distance between nodes%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Node_distance = zeros(No_nodes,No_nodes);
    for node = 1:No_nodes
        for j = 1:No_nodes
            Node_distance(node,j) = sqrt((abs(Nodes_position(node,1)-Nodes_position(j,1)))^2 + (abs(Nodes_position(node,2)-Nodes_position(j,2)))^2);
        end
    end
    
    % START OF INNER LOOP (Change of communication range)
    %---------------------------------------------------               
    for CommRange = 1:19
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %6. Step: Set communication range%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Max_comm_distance = Max_comm_d(1,CommRange); % Calculation of current comm range
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %7. Step: Flooding principles for specific communication range%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initialization of available nodes
        available_nodes = zeros(No_nodes,No_nodes);
        for i = 1:No_nodes
            j = i-1;
            while j ~= 0
                if Node_distance(i,j) <= Max_comm_distance
                    available_nodes(i,j) = 1;
                end
                j = j-1;
            end 
        end
        available_nodes = available_nodes + available_nodes.';
        next_source = 0;
        
        % Basic flooding %
        hops = 0;
        messages = 0;
        found = 0;
        source = Index_S_D(1,1);
        while found ~= 1
            if ~isempty(source)
               messages = messages + length(source);
            end
            hops = hops + 1;
            [~, next_source] = ind2sub(size(available_nodes(source,:)),find(available_nodes(source,:)));
            next_source = unique(next_source(:).');
            if any(next_source - Index_S_D(2,1) == 0) || hops == TTL
                found = 1;
            else
                source = next_source;
            end
        end
        Total_hops(CommRange,1) = Total_hops(CommRange,1) + hops;
        Total_messages(CommRange,1) = Total_messages(CommRange,1) + messages;
        
        % Improved flooding %
        hops = 0;
        messages = 0;
        found = 0;
        source = Index_S_D(1,1);
        is_unique = ones(No_nodes,1);
        while found ~= 1
            if ~isempty(source)
                source = setdiff(source,find(is_unique == 0));  % separate non-unique sources
                messages = messages + length(source);
                is_unique(source) = 0;                          % set current sources as non-unique
            end
            hops = hops + 1;
            [~, next_source] = ind2sub(size(available_nodes(source,:)),find(available_nodes(source,:)));
            next_source = unique(next_source);
            if any(next_source - Index_S_D(2,1) == 0) || hops == TTL
                found = 1;
            else
                source = next_source;
            end
        end
        Total_hops(CommRange,2) = Total_hops(CommRange,2) + hops;
        Total_messages(CommRange,2) = Total_messages(CommRange,2) + messages;
        
        % Gossiping %
        hops = 0;
        messages = 0;
        found = 0;
        source = Index_S_D(1,1);
        is_unique = ones(No_nodes,1);
        while found ~= 1
            if ~isempty(source)
                source = setdiff(source,find(is_unique == 0));              % improved flooding step
                source = nonzeros(source.*binornd(1,p,1,length(source)));   % gossiping step
                messages = messages + length(source);
                is_unique(source) = 0;                                      % improved flooding step
            end
            hops = hops + 1;
            [~, next_source] = ind2sub(size(available_nodes(source,:)),find(available_nodes(source,:)));
            next_source = unique(next_source(:).');
            if any(next_source - Index_S_D(2,1) == 0) || hops == TTL
                found = 1;
            else
                source = next_source;
            end
        end
        Total_hops(CommRange,3) = Total_hops(CommRange,3) + hops;
        Total_messages(CommRange,3) = Total_messages(CommRange,3) + messages;

    end
    
    % END of INNER LOOP
    %---------------------------------------------------
    disp(simulation_number)
    
end
% END OF OUTER LOOP
%---------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%8. Step: Average the results over K simulation drops%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Total_messages=Total_messages/Sim_drops;
Total_hops=Total_hops/Sim_drops;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%9. Step: Plot the results%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure (1)
hold on
grid on
axis([1 19 0 ceil(max(Total_messages(:,1)))]);

set(gca, 'fontweight', 'normal', ...
    'fontsize', 14)
xlabel('Maximum communication range [m]')
ylabel('Average number of messages')
a=plot(1:19,Total_messages(1:19,1),'b>-','MarkerSize',6,'LineWidth',1);
b=plot(1:19,Total_messages(1:19,2),'r*-','MarkerSize',6,'LineWidth',1);
c=plot(1:19,Total_messages(1:19,3),'g+-','MarkerSize',6,'LineWidth',1);
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19])
set(gca,'XTickLabel',{'10','15','20','25','30','35','40','45','50','55','60','65','70','75','80','85','90','95','100'})
legend([a,b,c],'Basic flooding','Improved flooding',['Gossiping (p=' num2str(p) ')'],'Location','northeast');
print('-f1', '-r300', '-dpng','Messages')

figure (2)
hold on
grid on
axis([1 19 0 ceil(max(Total_hops(:,1)))]);

set(gca, 'fontweight', 'normal', ...
    'fontsize', 14)
xlabel('Maximum communication range [m]')
ylabel('Average number of hops')
a=plot(1:19,Total_hops(1:19,1),'b>-','MarkerSize',6,'LineWidth',1);
b=plot(1:19,Total_hops(1:19,2),'r*-','MarkerSize',6,'LineWidth',1);
c=plot(1:19,Total_hops(1:19,3),'g+-','MarkerSize',6,'LineWidth',1);
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19])
set(gca,'XTickLabel',{'10','15','20','25','30','35','40','45','50','55','60','65','70','75','80','85','90','95','100'})
legend([a,b,c],'Basic flooding','Improved flooding',['Gossiping (p=' num2str(p) ')'],'Location','northeast');
print('-f2', '-r300', '-dpng','Hops')
hold off
