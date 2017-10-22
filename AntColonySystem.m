clear;clc;
% INITIALIZE DATA
%nodes = csvread('city10.csv')';
nodes = csvread('eil51.csv')';
%nodes = csvread('berlin52.csv')';

% PARAMETER SETTING
MAX_ITERATION = 2000;
alpha = 0.1;
beta = 2;
rho = 0.1;
q0 = 0.9;
number_of_ants = 10;

iteration = 1;
number_of_nodes = size(nodes, 2);
distance = dist(nodes);
visibility = 1./distance;

shortest_distance_each_iteration = zeros(MAX_ITERATION, 1);
shortest_distance_all_iteration = zeros(MAX_ITERATION, 1);

visited_nodes = zeros(1, number_of_nodes);
unvisited_nodes = colon(1, number_of_nodes);

% Compute Lnn
Lnn = 0;
for n = 1:number_of_nodes
    if n==1
        visited_nodes(n) = floor(rand()*number_of_nodes+1);
        tempLnn = visited_nodes;
        tempLnn(tempLnn == 0) = [];
        unvisited_nodes(unvisited_nodes == tempLnn) = [];
    else
        distance_temp = zeros(1, size(unvisited_nodes, 2));
        for i = 1:size(unvisited_nodes, 2)
            if visited_nodes(n-1) ~= unvisited_nodes(i)
                distance_temp(unvisited_nodes(i)) = distance(visited_nodes(n-1), unvisited_nodes(i));
            end
        end
        distance_temp(distance_temp == 0) = NaN;
        founded = find(distance_temp == min(distance_temp));
        Lnn = Lnn + min(distance_temp);
        visited_nodes(n) = founded(1);
        unvisited_nodes(find(unvisited_nodes == founded(1))) = [];
    end
end % end of Compute Lnn

initial_tao = 1/(number_of_nodes*Lnn);
tao = ones(number_of_nodes, number_of_nodes) * initial_tao;

top_global_ant_path = zeros(1, number_of_nodes+1);
top_global_ant_distance = 0;

top_local_ant_path = zeros(1, number_of_nodes+1);
top_local_ant_distance = 0;

while iteration <= MAX_ITERATION
    s=1; % s: tabu list index -> s-th node will be filled
    tabu = zeros(number_of_ants, number_of_nodes+1); % tabu for saving the route of each ant
    tabu(:, 1) = floor((rand(number_of_ants, 1)*number_of_nodes)+1); % random initial node of ant

    % fill tabu for the next node
    while s <= number_of_nodes
        s = s + 1; % s = tabu index to be filled
        for k = 1:number_of_ants % looping according to total ants | canculate the probability of each ant
            if s == (number_of_nodes+1) % last tabu, fill it with initial node
                tabu(:, end) = tabu(:, 1);
            else
                % determine the next node to be visited
                q = rand();
                probability = zeros(1, number_of_nodes);
                temp = colon(1, number_of_nodes); % matrix for unvisited node
                tabu1 = tabu(k, 1:end-1); % create clone of tabu
                tabu1(tabu1 == 0) = []; % remove visited node
                temp(tabu1) = []; % remove visited node
                
                for x = 1:size(temp, 2)
                    probability(1, temp(x)) = tao(tabu(k, s-1), temp(x)) * visibility(tabu(k, s-1), temp(x))^beta;
                end
                
                if q <= q0 % exploitation, when q <= q0
                    founded = find(probability == max(probability));
                    tabu(k, s) = founded(floor((rand()*size(founded, 2))+1));
                else % explorasi, when q > q0
                    probability = probability./sum(probability);
                    % choose the next node based on probability
                    sorted = sort(probability);
                    range = cumsum(sorted);
                    randomResult = rand();
                    founded1 = find(probability == sorted(max(find(range <= randomResult))+1));
                    tabu(k, s) = founded1(floor((rand()*size(founded1, 2))+1)); % random when there're the same probabilities                    
                end
            end % k-th ant and has been finished moving to the city
            
            % LOCAL UPDATE -> tao = ((1-rho) .* tao) + (rho .* initial_tao);
            tao(tabu(k, s-1), tabu(k, s)) = ((1-rho) .* tao(tabu(k, s-1), tabu(k, s))) + (rho .* initial_tao);
            tao(tabu(k, s), tabu(k, s-1)) = tao(tabu(k, s-1), tabu(k, s));
        end % all ants' s+1 tabu has been filled
    end
    
    % calculate the distance
    distance_of_ants_tour = zeros(number_of_ants, 1);
    for k=1:number_of_ants
        for x = 1:number_of_nodes
            distance_of_ants_tour(k, 1) = distance_of_ants_tour(k, 1) + round(distance(tabu(k, x), tabu(k, x+1)));
        end
    end
    
    % update 'shortest_path' with the shortest_path of current iteration
    founded = find(distance_of_ants_tour == min(distance_of_ants_tour));
    if iteration == 1 || (min(distance_of_ants_tour) < shortest_distance_all_iteration(iteration-1))
        shortest_distance_all_iteration(iteration) = min(distance_of_ants_tour);
        
        % update top global ant
        top_global_ant_path = tabu(founded(1), :);
        top_global_ant_distance = min(distance_of_ants_tour);
    else
        shortest_distance_all_iteration(iteration) = shortest_distance_all_iteration(iteration-1);
    end
    
    % update top local ant
    top_local_ant_path = tabu(founded(1), :);
    top_local_ant_distance = min(distance_of_ants_tour);
    
    % update shortest distance of each iteration
    shortest_distance_each_iteration(iteration) = min(distance_of_ants_tour);
    
    % GLOBAL UPDATE pheromone from top ant  
    for x = 1:number_of_nodes
        % LOCAL BEST ANT -> tao = ((1-alpha).*tao) + (alpha .* delta_tao);
         tao(top_local_ant_path(x), top_local_ant_path(x+1)) = ((1-alpha).*tao(top_local_ant_path(x), top_local_ant_path(x+1))) + (alpha .* (1/top_local_ant_distance));
         tao(top_local_ant_path(x+1), top_local_ant_path(x)) = tao(top_local_ant_path(x), top_local_ant_path(x+1));        
    end
    fprintf('iteration: %d/%d | shortest distance of all iteration: %d \n', iteration, MAX_ITERATION, shortest_distance_all_iteration(iteration));
    iteration = iteration + 1;
end

close
f = figure;
subplot(1, 2, 1);
plot([1:MAX_ITERATION]', shortest_distance_each_iteration(:, 1)');
title('Shortest Distance Each Iteration');
xlabel('Iterations');
ylabel('Shortest Distance');

subplot(1, 2, 2);
plot([1:MAX_ITERATION]', shortest_distance_all_iteration(:, 1)');
title('Shortest Distance All Iteration');
xlabel('Iterations');
ylabel('Shortest Distance');
drawnow
set(f, 'position', [200, 200, 700, 300]);
disp(top_global_ant_path);
disp(shortest_distance_all_iteration(MAX_ITERATION));