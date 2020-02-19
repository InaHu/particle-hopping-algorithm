% This program solves the time evolution of both vesicle transport and
% pool concentration with a Particle-Hopping-Ansatz. The algorithm is
% explained in more detail in the paper 
% "On the Role of Vesicle Transport in Neurite Growth: Modelling and Experiments"
% Please read the README-file that contains important information on how 
% to handle this code.

%% Load initial conditions
close all   
clear variables

% Initialise all neurites as structure arrays
NumberOfNeurites = 2;
N1 = struct;     

% Typical values (all parameters that have a tilde in the paper, see Chapter 3.1 'Scaling')
N1.typicalTime = 100;                   
N1.typicalLength = 50;
N1.typicalDiffusion = 10^(-1); 
N1.typicalInflux = 1;
N1.typicalOutflux = 10^(-1);
N1.typicalPotential = 1;
N1.typicalConcentration = 15;

% Calculate Scaling Parameters (see Chapter 3.1 'Scaling')
lambda_eps = N1.typicalTime*N1.typicalDiffusion/(2*N1.typicalLength^2);
lambda_in = (N1.typicalTime*N1.typicalInflux)/(2*N1.typicalLength*N1.typicalConcentration);       
lambda_out = (N1.typicalTime*N1.typicalOutflux)/(2*N1.typicalLength);

N1.lambda_in = lambda_in;
N1.lambda_out = lambda_out;
N1.lambda_eps = lambda_eps;

% Second neurite hat the same scaling parameters
N2 = N1;

% Space discretisation
n = 400;
N1.L = 1;                                  % Length of Neurite 1
N2.L = 0.3;
%N2.L = 0.9;
N1.h = 1/n;                                % Space step size                                
N2.h = 1/n;
N1.x = linspace(0,N1.L,n*N1.L)';           % Domain of Neurite 1            
N2.x = linspace(0,N2.L,n*N2.L)';            

% Time discretisation
T = 15;                                   % End time
m = T*10000;
tau = T/m;                                 % Time step size
N1.tau = tau;   
N2.tau = tau;                            
tt = linspace(0,T,ceil(T/N1.tau))';  

% Load initial pool Cconcentrations and define their maximum sizes (see Chapter 5.1);
N1.Lambda_som = 0.12;                        % Total Mass of Vesicles in the Soma
N2.Lambda_som = N1.Lambda_som;
N1.Lambda_tip = 0.0015;                      % Total Mass of Vesicles in the growth cone of N1
N2.Lambda_tip = 0.0015;                      % Total Mass of Vesicles in the growth cone of N2

N1.Lambda_som_max = 0.175;
N2.Lambda_som_max = N1.Lambda_som_max;
N1.Lambda_tip_max = 0.0029;               
N2.Lambda_tip_max = 0.0029;

SaveLambda_som = N1.Lambda_som; 
SaveN1Lambda_tip = N1.Lambda_tip; 
SaveN2Lambda_tip = N2.Lambda_tip;

% Initial influx- and outflux values in all neurites (see Chapter 5.1)
no_noflux_1 = 1;                          % Set to 0 to have no flux in neurite 1       
no_noflux_2 = 1;                          % Set to 0 to have no flux in neurite 2
N1.alpha_a = no_noflux_1*0.4;               
N1.alpha_r = no_noflux_1*0.2;         
N2.alpha_a = no_noflux_2*0.4;
N2.alpha_r = no_noflux_2*0.2;

N1.beta_a = no_noflux_1*15;
N1.beta_r = no_noflux_1*15;
N2.beta_a = no_noflux_2*15;
N2.beta_r = no_noflux_2*15;

% Further parameters (diffusion constants and potentials) (see Chapter 5.1)
N1.eps_a = 0.01;
N2.eps_a = N1.eps_a;
N1.eps_r = 0.01;
N2.eps_r = N1.eps_r;

N1.V_a = 1.75.*N1.x;           
N1.V_r = -1.5.*N1.x;
N2.V_a = 1.75.*N2.x;
N2.V_r = -1.5.*N2.x;

% Initial concentration of a = ANT and r = RET (see Chapter 5.1)
N1.a0 = 0.*N1.x + 0.1;
N1.r0 = 0.*N1.x + 0.1;
N2.a0 = 0.*N2.x + 0.1;
N2.r0 = 0.*N2.x + 0.1;

% Initialise vectors that save development of different numbers
Development_Lambda_som = zeros(m,1);
Development_Lambda_tipN1 = zeros(m,1);
Development_Lambda_tipN2 = zeros(m,1);
Development_MassWholeSystem = zeros(m,1);     

% Define Colors
redN = [.8 .2941  .0862];
greenN = [ .1647 .6313 .5960];
greyN = [ 0.4 0.4 0.4];

%% Check initial concentrations (non-negativity and density constraint)

N1.p0 = N1.a0 + N1.r0;                  % total concentration in neurite 1
N2.p0 = N2.a0 + N2.r0;

if sum(N1.p0 > 1) + sum(N2.p0 > 1) > 0
    error('The chosen initial value of the total vesicle concentration in one neurite exceeds the density constraint.')
end
if sum(N1.a0 < 0) + sum(N1.r0 < 0) + sum(N2.a0 < 0) + sum(N2.r0 < 0) > 0
    error('One of the initial vesicle concentrations in a neurite is negative.')
end
if N1.Lambda_som > N1.Lambda_som_max || N1.Lambda_tip > N1.Lambda_tip_max || N2.Lambda_tip > N2.Lambda_tip_max
    error('One of the initial pool concentration exceeds the maximum concentration.')
end

% Calculate initial mass 
TotalMass = N1.Lambda_som  + N1.Lambda_tip  +  N2.Lambda_tip  ... 
         + (sum(N1.r0(:)) + sum(N1.a0(:)) + sum(N2.r0(:)) + sum(N2.a0(:)))*N2.h;
disp(['Total mass at the beginning is = ' num2str(TotalMass), '.']);

%% Plot Initial Concentration in Both Neurites, Neurite 1 is plotted reverse
fig = figure('Position', [200, 100, 1200, 600]);     
  
plotlefttip = subplot(1,5,1);
hl = bar(N1.Lambda_tip, 'FaceColor', greyN);
ylim([0 N1.Lambda_tip_max ]);
plotlefttip.PlotBoxAspectRatio(2) = 2.5;
plotlefttip.Position(1) = 0.02;
plotlefttip.Position(2) = .5;
plotlefttip.Position(4) = 0.4;
set(gca,'FontSize',15,'FontWeight','bold')

fs1 = subplot(1,5,2);
h1 = plot(N1.x, N1.a0, N1.x, N1.r0, N1.x, 1-N1.a0-N1.r0 ,'k--','Linewidth',2);
h1(1).Color = greenN;
h1(2).Color = redN;
fs1.Position = [0.1000    0.5000    0.3237    0.4150];

set(fs1, 'Xdir', 'reverse')
title('Neurite 1')
axis([0 N1.L 0 1]),
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultAxesFontName', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
legend('ANT','RET', '1-$\rho$');
xlabel('Points on the domain $\Omega_1$');
ylabel('Concentration');
axis square
set(gca,'FontSize',15,'FontWeight','bold')

plotsoma = subplot(1,5,3);
hs = bar(N1.Lambda_som, 'FaceColor', greyN);
title('$\Lambda_{som}$')
ylim([0 N1.Lambda_som_max]);
plotsoma.PlotBoxAspectRatio(2) = 2.5;
plotsoma.Position(1) = 0.36;
plotsoma.Position(2) = 0.5;
plotsoma.Position(4) = 0.4;
set(gca,'FontSize',15,'FontWeight','bold')

fs2 = subplot(1,5,4);
h2 = plot(N2.x, N2.a0, N2.x, N2.r0, N2.x, 1-N2.a0-N2.r0 ,'k--','Linewidth',2);
h2(1).Color = greenN;
h2(2).Color = redN;
fs2.Position = [0.4400    0.5000    0.3237    0.4150];
axis([0 N2.L 0 1]);
title('Neurite 2')
axis square
legend('ANT','RET', '1-$\rho$');
xlabel('Points on the domain $\Omega_2$');
ylabel('Concentration');
set(gca,'FontSize',15,'FontWeight','bold','Fontangle','Oblique')

plotrighttip = subplot(1,5,5);
hr = bar(N2.Lambda_tip, 'FaceColor', greyN);
ylim([0 N2.Lambda_tip_max ]);
plotrighttip.PlotBoxAspectRatio(2) = 2.5;
plotrighttip.Position(1) = 0.7;
plotrighttip.Position(2) = 0.5;
plotrighttip.Position(4) = 0.4;
set(gca,'FontSize',15,'FontWeight','bold')

drawnow;
exportFigure(['Images/TimeEvolution2B',num2str(0),'.pdf'], fig);          % command provided by helperFiles package, see README 
disp('The figure of the initial datum was printed.');

tic

%% Solve the equation
z = 1;                                                                     % Counter for grafic enumeration
for i = 1:m
    % Update Concentration in Neurons
    N1old = N1;
    [N1a1, N1r1] = UpdateConcentrationParticleHopping(N1);    
    N1.a0 = N1a1;
    N1.r0 = N1r1;
    N1p1 = N1a1 + N1r1;
	
    N2old = N2;
    [N2a1, N2r1] = UpdateConcentrationParticleHopping(N2);    
    N2.a0 = N2a1;
    N2.r0 = N2r1;
    N2p1 = N2a1 + N2r1;
    
    % Plot Updaten
    if (mod(i,1000)==0)

        h1(1).YData = N1a1;    
        h1(2).YData = N1r1;
        h1(3).YData = 1-N1a1-N1r1;

        h2(1).YData = N2a1;          
        h2(2).YData = N2r1;
        h2(3).YData = 1-N2a1-N2r1;
        
        hs.YData = N1.Lambda_som;
        hl.YData = N1.Lambda_tip;
        hr.YData = N2.Lambda_tip;
       
        drawnow
        
        disp(['Iteration: ', num2str(i)]);
        if (mod(i,10000)==0)
            %exportFigure(['Images/TimeEvolution2B',num2str(z),'.pdf'], fig);         % command provided by helperFiles package
            %disp(['A figure was printed at t = ',num2str(z), '.']);
            z = z + 1;
        end
    end
    
    % Save Density in Pools
    Development_Lambda_som(i) = N1.Lambda_som;
    Development_Lambda_tipN1(i) = N1.Lambda_tip;
    Development_Lambda_tipN2(i) = N2.Lambda_tip;
    
    % Update Concentration in Pools (see Chapter 5)
    N1.Lambda_som = N1.Lambda_som ...                               % Old concentration
        + N1.tau*(lambda_out*N1.beta_r*(1 - N1.Lambda_som/N1.Lambda_som_max)*(N1old.r0(1)) ...      % Retrograte entering from Neuron 1
            + lambda_out*N2.beta_r*(1 - N2.Lambda_som/N2.Lambda_som_max)*(N2old.r0(1)) ...               % Retrograte entering from Neuron 2
            - lambda_in*N1.alpha_a*(N1old.r0(1)-N1old.a0(1))*(N1.Lambda_som/N1.Lambda_som_max)*(1-N1old.r0(1)-N1old.a0(1)) ...                          % Anterograde existing pool into Neuron 1 
            - lambda_in*N2.alpha_a*(N2old.r0(1)-N2old.a0(1))*(N2.Lambda_som/N2.Lambda_som_max)*(1-N2old.r0(1)-N2old.a0(1)));                     % Anterograde existing pool into Neuron 2
            
    N2.Lambda_som = N1.Lambda_som;
    N1.Lambda_tip = N1.Lambda_tip ...                                           % Tip of Neuron 1
        + N1.tau*(- lambda_in*N1.alpha_r*(N1.Lambda_tip/ N1.Lambda_tip_max)*(1-N1old.r0(end)-N1old.a0(end)) ...                 % Vesicles exiting pool Neuron 1
            + lambda_out*N1.beta_a*N1old.a0(end)*(1 - N1.Lambda_tip/N1.Lambda_tip_max )); % Vesicles entering pool Neuron 1
    N2.Lambda_tip = N2.Lambda_tip ...                                           % Tip of Neuron 2
        + N2.tau*(- lambda_in*N2.alpha_r*(N2.Lambda_tip/N2.Lambda_tip_max)*(1-N2old.r0(end)-N2old.a0(end)) ...                 % Vesicles exiting pool Neuron 2
            + lambda_out*N2.beta_a*N2old.a0(end)*(1 - N2.Lambda_tip/N2.Lambda_tip_max));  

    % Calculate total mass (the system is mass conserving, so the mass will not change)
    Development_MassWholeSystem(i) = N1.Lambda_som + N1.Lambda_tip + N2.Lambda_tip  ... 
        + N2.h*(sum(N1r1(:)) + sum(N1a1(:)) + sum(N2r1(:)) + sum(N2a1(:))); 
end

%% Plot Pool Density Evolution in the pools
figure(5);
subplot(1,3,1);
Pool1 = plot(tt, Development_Lambda_som);
Pool1.Color = greyN;
axis square
title('$\Lambda_{som}$');
xlabel('Time $t$');
xlim([0 T]);
ylabel('Concentration');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultAxesFontName', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(gca,'FontSize',12,'FontWeight','bold');

subplot(1,3,2);
Pool2 = plot(tt, Development_Lambda_tipN1);
Pool2.Color = greyN;
axis square
title('$\Lambda_{N1}$');
xlabel('Time $t$');
ylabel('Concentration');
xlim([0 T]);
ylim([0 0.002]);
set(gca,'FontSize',12,'FontWeight','bold');

subplot(1,3,3),
Pool3 = plot(tt, Development_Lambda_tipN2);
Pool3.Color = greyN;
axis square
title('$\Lambda_{N2}$');
xlabel('Time $t$');
ylabel('Concentration');
xlim([0 T]);
ylim([0 0.002]);
set(gca,'FontSize',12,'FontWeight','bold');

resizeFigure(gcf, [270,850]);
exportFigure('Images/DevelopmentPoolConcentration2B.pdf', gcf);                      % command provided by helperFiles package
disp('A figure of the Development of the Concentration in the pools was printed.')

toc

