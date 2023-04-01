clc, clear, close all;

%% Problem Definition

model=CreateModel();

model.n=5;  % number of Handle Points

CostFunction=@(x) MyCost(x,model); 

nVar=model.n;

VarSize=[1 nVar];       % Decision Variables Matrix Size

VarMin.x=model.xmin;           
VarMin.y=model.ymin;           % Lower Bound of Variables
VarMax.x=model.xmax;           % Upper Bound of Variables
VarMax.y=model.ymax;         
 
%% Firefly Algorithm Parameters

MaxIt=200;         % Maximum Number of Iterations
nPop=60;            % Number of Fireflies (Swarm Size)
gamma=1;            % Light Absorption Coefficient
beta0=2;            % Attraction Coefficient Base Value
alpha=0.2;          % Mutation Coefficient
alpha_damp=0.98;    % Mutation Coefficient Damping Ratio

delta.x = 0.05*(VarMax.x-VarMin.x); 
delta.y = 0.05*(VarMax.y-VarMin.y); 

m=2;

%%% x %%%
if isscalar(VarMin.x) && isscalar(VarMax.x)
    dmax.x = (VarMax.x-VarMin.x)*sqrt(nVar);
else
    dmax.x = norm(VarMax.x-VarMin.x);
end

%%% y %%%
if isscalar(VarMin.y) && isscalar(VarMax.y)
    dmax.y = (VarMax.y-VarMin.y)*sqrt(nVar);
else
    dmax.y = norm(VarMax.y-VarMin.y);
end

%% Initialization

% Empty Firefly Structure
firefly.Position=[];
firefly.Cost=[];
firefly.Sol=[];

% Initialize Population Array
pop = repmat(firefly,nPop,1);

% Initialize Best Solution Ever Found
BestSol.Cost=inf;

% Create Initial Fireflies
for i=1:nPop
    if i > 1
        pop(i).Position = CreateRandomSolution(model);
    else
        % Straight line from source to destination
        xx = linspace(model.xs, model.xt, model.n+2);
        yy = linspace(model.ys, model.yt, model.n+2);
        pop(i).Position.x = xx(2:end-1);
        pop(i).Position.y = yy(2:end-1);
   end
    
    [pop(i).Cost, pop(i).Sol] = CostFunction(pop(i).Position);
    
    pop(i).Best.Position=pop(i).Position;
    pop(i).Best.Cost=pop(i).Cost;
    pop(i).Best.Sol=pop(i).Sol;
    
    if pop(i).Best.Cost <= BestSol.Cost
       BestSol=pop(i).Best;
    end
   
end

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% Firefly Algorithm Main Loop

for it=1:MaxIt
    
    newpop=repmat(firefly,nPop,1);
    for i=1:nPop
        newpop(i).Best.Cost= inf;
        
        for j=1:nPop
            if pop(j).Cost < pop(i).Cost
                rij = norm(pop(i).Position.x - pop(j).Position.x)/dmax.x;
                beta = beta0 * exp(-gamma*rij^m);
                e = delta.x * randn(VarSize);
                
                newsol.Position.x = pop(i).Position.x ...
                                + beta*rand(VarSize).*(pop(j).Position.x-pop(i).Position.x) ...
                                + alpha*e;
                
                newsol.Position.x = max(newsol.Position.x,VarMin.x);
                newsol.Position.x = min(newsol.Position.x,VarMax.x);
                
              %%% y %%%  
              
                rij = norm(pop(i).Position.y - pop(j).Position.y)/dmax.y;
                beta = beta0 * exp(-gamma*rij^m);
                e = delta.y * randn(VarSize);
                newsol.Position.y = pop(i).Position.y ...
                                + beta * rand(VarSize).*(pop(j).Position.y - pop(i).Position.y) ...
                                + alpha * e;
                
                newsol.Position.y = max(newsol.Position.y,VarMin.y);
                newsol.Position.y = min(newsol.Position.y,VarMax.y);
                
                [newsol.Cost, newsol.Sol] = CostFunction(newsol.Position);
                newsol.Best.Position=newsol.Position;
                newsol.Best.Cost=newsol.Cost;
                newsol.Best.Sol=newsol.Sol;
                
                if newsol.Best.Cost <= newpop(i).Best.Cost
                    newpop(i).Best = newsol.Best;
                    if newpop(i).Best.Cost <= BestSol.Cost
                        BestSol = newpop(i).Best;
                    end
                end
                
            end
        end
    end
    
    % Merge
    pop=[pop 
        newpop];  %#ok
    
    % Sort
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    
    % Truncate
    pop=pop(1:nPop);
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Damp Mutation Coefficient
    alpha = alpha*alpha_damp;
    
    if BestSol.Sol.IsFeasible
        Flag='*';
    else
        Flag=[', Violation = ' num2str(BestSol.Sol.Violation)];
    end
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    % Plot Solution
    figure(1);
    PlotSolution(BestSol.Sol, model);
    pause(0.01); 
end

%% Results

figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
