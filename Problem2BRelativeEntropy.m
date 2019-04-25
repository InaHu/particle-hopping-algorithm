% This program solves the time evolution of problem 2B with a particle hopping ansatz
%% Load initial conditions
close all   
clear variables

% Initialise all neurites as structure arrays
NumberOfNeurites = 2;
N1 = struct;     

% Initialise parameters for scaling
N1.typicalTime = 100;
N1.typicalLength = 50;
N1.typicalDiffusion = 10^(-1); 
N1.typicalInOutFlux = 0.1;
N1.typicalPotential = 1;
N1.typicalConcentration = 15;

N2 = N1;

% Space discretisation
n = 400;
N1.L = 1;
N2.L = 0.3;
N1.dx = 1/n;                                % or h
N2.dx = 1/n;
N1.x = linspace(0,N1.L,n*N1.L)';            
N2.x = linspace(0,N2.L,n*N2.L)';            

% Time discretisation
T = 15;                                     % End time
m = T*10000;
dt = T/m;
N1.dt = dt;   
N2.dt = dt;                            
tt = linspace(0,T,ceil(T/N1.dt))';  

% Load initial pool Cconcentrations and define their maximum sizes;
N1.Lambda_som = 1.4;                        % Total Mass of Vesicles in the Soma
N2.Lambda_som = N1.Lambda_som;
N1.Lambda_tip = 0;                          % Total Mass of Vesicles in the growth cone of N1
N2.Lambda_tip = 0;                          % Total Mass of Vesicles in the growth cone of N2

N1.Lambda_som_max = 1.4;
N2.Lambda_som_max = N1.Lambda_som_max;
N1.Lambda_tip_max = 0.02;               
N2.Lambda_tip_max = 0.02;

SaveLambda_som = N1.Lambda_som; 
SaveN1Lambda_tip = N1.Lambda_tip; 
SaveN2Lambda_tip = N2.Lambda_tip;

% Initial influx- and outflux values in all neurites
no_noflux_1 = 1;                          % Set to 0 to have no flux in neurite 1
no_noflux_2 = 1;                          % Set to 0 to have no flux in neurite 2
N1.alpha_a = no_noflux_1*8;               % 5 means that 50% would enter 
N1.alpha_r = no_noflux_1*5.5;         
N2.alpha_a = no_noflux_2*8;
N2.alpha_r = no_noflux_2*5.5;

N1.beta_a = no_noflux_1*8;
N1.beta_r = no_noflux_1*5.5;
N2.beta_a = no_noflux_2*8;
N2.beta_r = no_noflux_2*5.5;

% Further parameters (diffusion constants and potentials)
N1.eps_a = 0.05;
N2.eps_a = N1.eps_a;
N1.eps_r = 0.05;
N2.eps_r = N1.eps_r;

N1.V_a = 1.75.*N1.x;           
N1.V_r = -1.5.*N1.x;
N2.V_a = 1.75.*N2.x;
N2.V_r = -1.5.*N2.x;

% With eps we can calcultate the Scaling for the Pool
K = (N1.typicalInOutFlux)/(2*N1.dx*N1.typicalConcentration);

% Initial concentration of a = ANT and r = RET
N1.a0 = 0.*N1.x + 0.2;
N1.r0 = 0.*N1.x + 0.21;
N2.a0 = 0.*N2.x + 0.2;
N2.r0 = 0.*N2.x + 0.21;

SaveN1a0 = N1.a0; 
SaveN1r0 = N1.r0; 
SaveN2a0 = N2.a0;
SaveN2r0 = N2.r0;

% Initialise vectors that save development of different numbers
Development_Lambda_som = zeros(m,1);
Development_Lambda_tipN1 = zeros(m,1);
Development_Lambda_tipN2 = zeros(m,1);
Development_MassWholeSystem = zeros(m,1);
Development_Entropy = zeros(m,1);

%% Check initial concentrations

N1.p0 = N1.a0 + N1.r0;
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

% % Calculate initial mass ACHTUNG AUFPASSEN WIE LAMBDA DEFINIERT
% TotalMass = N1.Lambda_som  + N1.Lambda_tip  +  N2.Lambda_tip  ... 
%          + (sum(N1.r0(:)) + sum(N1.a0(:)) + sum(N2.r0(:)) + sum(N2.a0(:)))*N2.dx;
% disp(['Total mass at the beginning is = ' num2str(TotalMass), '.']);

%% Plot Initial Concentration in Both Neurites, Neurite 1 is plotted reverse
%fig = figure('Position', [0, 0, 200, 200]);        % for Jans Computer
fig = figure('Position', [200, 100, 1200, 600]);     % for Inas Computer
  
plotlefttip = subplot(1,5,1);
hl = bar(N1.Lambda_tip);
title('$\Lambda_{N1}$')
ylim([0 N1.Lambda_tip_max ]);
plotlefttip.PlotBoxAspectRatio(2) = 2.5;
plotlefttip.Position(1) = 0.02;
plotlefttip.Position(2) = .5;
plotlefttip.Position(4) = 0.4;
set(gca,'FontSize',15,'FontWeight','bold')

fs1 = subplot(1,5,2);
h1 = plot(N1.x, N1.a0,'blue', N1.x, N1.r0,'red', N1.x, 1-N1.a0-N1.r0 ,'k--','Linewidth',2);
% fs1.Position = [0.18000    0.5000    0.3237    0.4150];
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
hs = bar(N1.Lambda_som);
title('$\Lambda_{som}$')
ylim([0 N1.Lambda_som_max]);
plotsoma.PlotBoxAspectRatio(2) = 2.5;
% plotsoma.Position(1) = 0.5;
plotsoma.Position(1) = 0.36;
plotsoma.Position(2) = 0.5;
plotsoma.Position(4) = 0.4;
set(gca,'FontSize',15,'FontWeight','bold')

fs2 = subplot(1,5,4);
h2 = plot(N2.x, N2.a0, 'blue', N2.x, N2.r0,'red',N2.x, 1-N2.a0-N2.r0 ,'k--','Linewidth',2);
fs2.Position = [0.4400    0.5000    0.3237    0.4150];
axis([0 N2.L 0 1]);
title('Neurite 2')
axis square
legend('ANT','RET', '1-$\rho$');
xlabel('Points on the domain $\Omega_2$');
ylabel('Concentration');
set(gca,'FontSize',15,'FontWeight','bold','Fontangle','Oblique')

plotrighttip = subplot(1,5,5);
hr = bar(N2.Lambda_tip);
title('$\Lambda_{N2}$')
ylim([0 N2.Lambda_tip_max ]);
plotrighttip.PlotBoxAspectRatio(2) = 2.5;
plotrighttip.Position(1) = 0.7;
plotrighttip.Position(2) = 0.5;
plotrighttip.Position(4) = 0.4;
set(gca,'FontSize',15,'FontWeight','bold')

fs1.Position = [0.1000    0.5000    0.3237    0.4150];

drawnow;

%exportFigure(['TimeEvolution2B',num2str(0),'.pdf'], fig);
disp('The figure of the initial datum was printed.');

tic;

% %% Calculate Stationary Solution
% 
% disp('Wait until program has calculated the stationary solutions.');
% 
% for i = 1:(m*2)
%     % Update Concentration in Neurons
%     [N1a1, N1r1, ~] = UpdateConcentrationParticleHopping(N1);
%     N1old = N1;
%     N1.a0 = N1a1;
%     N1.r0 = N1r1;
%     N1p1 = N1a1 + N1r1;
% 	
%     [N2a1, N2r1, ~] = UpdateConcentrationParticleHopping(N2);    
%     N2old = N2;
%     N2.a0 = N2a1;
%     N2.r0 = N2r1;
%     N2p1 = N2a1 + N2r1;
%       
%     % Update Concentration in Pools 
%     N1.Lambda_som = N1.Lambda_som ...                               % Old concentration
%         + N1.dt*K*(N1.beta_r*(1 - N1.Lambda_som/N1.Lambda_som_max)*(N1r1(1)) ...      % Retrograte entering from Neuron 1
%             + N2.beta_r*(1 - N2.Lambda_som/N2.Lambda_som_max)*(N2r1(1)) ...               % Retrograte entering from Neuron 2
%             - N1.alpha_a*(N1.Lambda_som/N1.Lambda_som_max)*(1-N1r1(1)-N1a1(1)) ...                          % Anterograde existing pool into Neuron 1 
%             - N2.alpha_a*(N2.Lambda_som/N2.Lambda_som_max)*(1-N2r1(1)-N2a1(1))); ...                        % Anterograde existing pool into Neuron 2
%             
%     N2.Lambda_som = N1.Lambda_som;
%     N1.Lambda_tip = N1.Lambda_tip ...                                           % Tip of Neuron 1
%         + N1.dt*K*(- N1.alpha_r*(N1.Lambda_tip/ N1.Lambda_tip_max)*(1-N1r1(end)-N1a1(end)) ...                 % Vesicles exiting pool Neuron 1
%             + N1.beta_a*N1a1(end)*(1 - N1.Lambda_tip/N1.Lambda_tip_max ));                % Vesicles entering pool Neuron 1
%     N2.Lambda_tip = N2.Lambda_tip ...                                           % Tip of Neuron 2
%         + N1.dt*K*(- N2.alpha_r*(N2.Lambda_tip/N1.Lambda_tip_max)*(1-N2r1(end)-N2a1(end)) ...                 % Vesicles exiting pool Neuron 2
%             + N2.beta_a*N2a1(end)*(1 - N2.Lambda_tip/N2.Lambda_tip_max));  
%         
%     %disp(['sum(N1r1(:) = '  num2str(sum(N1r1(:))) ', sum(N1a1(:)) = ' num2str(sum(N1a1(:))) ', Lambda_som = ' num2str(N1.Lambda_som) ', N1.Lambda_tip = ' num2str(N1.Lambda_tip)]);
%     %disp(['What leaves soma = ' num2str(- N1.Scaling_Pools*N1.dt*N1.alpha_a*N1.Lambda_som*(1-N1old.r0(1)-N1old.a0(1))) ', (1-\rho)= ' num2str(1-N1old.r0(1)-N1old.a0(1))]);
% end
% 
% % Save stationary solutions
% N1.a_infinity = N1a1;
% N1.r_infinity = N1r1;
% N2.a_infinity = N2a1;
% N2.r_infinity = N2r1;
% 
% N1.Lambda_som_infinity = N1.Lambda_som;
% N2.Lambda_som_infinity = N1.Lambda_som_infinity;
% N1.Lambda_tip_infinity = N1.Lambda_tip;
% N2.Lambda_tip_infinity = N2.Lambda_tip;
% 
% 
% %% Reset Initial Concentrations
% 
% N1.Lambda_som = SaveLambda_som;                     
% N2.Lambda_som = SaveLambda_som;
% N1.Lambda_tip = SaveN1Lambda_tip;                   
% N2.Lambda_tip = SaveN2Lambda_tip;    
% 
% N1.a0 = SaveN1a0;
% N1.r0 = SaveN1r0;
% N2.a0 = SaveN2a0;
% N2.r0 = SaveN2r0;

%% Solve the equation

for i = 1:m
    % Update Concentration in Neurons
    [N1a1, N1r1, ~] = UpdateConcentrationParticleHopping(N1);    
    N1.a0 = N1a1;
    N1.r0 = N1r1;
    N1p1 = N1a1 + N1r1;
	
    [N2a1, N2r1, ~] = UpdateConcentrationParticleHopping(N2);    
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
        if (mod(i,20000)==0)
            % print(['Images/TimeEvolution2B',num2str(i*dt*10)], '-dpdf', '-bestfit');
            % exportFigure(['TimeEvolution2B',num2str(i*dt*10),'.pdf'], fig);
            disp(['A figure was printed at t = ',num2str(i*dt*10),num2str(N2.L), '.']);
        end
    end
    
    % Save Density in Pools
    Development_Lambda_som(i) = N1.Lambda_som;
    Development_Lambda_tipN1(i) = N1.Lambda_tip;
    Development_Lambda_tipN2(i) = N2.Lambda_tip;
    
    % Update Concentration in Pools 
    N1.Lambda_som = N1.Lambda_som ...                               % Old concentration
        + N1.dt*K*(N1.beta_r*(1 - N1.Lambda_som/N1.Lambda_som_max)*(N1r1(1)) ...      % Retrograte entering from Neuron 1
            + N2.beta_r*(1 - N2.Lambda_som/N2.Lambda_som_max)*(N2r1(1)) ...               % Retrograte entering from Neuron 2
            - N1.alpha_a*(N1.Lambda_som/N1.Lambda_som_max)*(1-N1r1(1)-N1a1(1)) ...                          % Anterograde existing pool into Neuron 1 
            - N2.alpha_a*(N2.Lambda_som/N2.Lambda_som_max)*(1-N2r1(1)-N2a1(1))); ...                        % Anterograde existing pool into Neuron 2
            
    N2.Lambda_som = N1.Lambda_som;
    N1.Lambda_tip = N1.Lambda_tip ...                                           % Tip of Neuron 1
        + N1.dt*K*(- N1.alpha_r*(N1.Lambda_tip/ N1.Lambda_tip_max)*(1-N1r1(end)-N1a1(end)) ...                 % Vesicles exiting pool Neuron 1
            + N1.beta_a*N1a1(end)*(1 - N1.Lambda_tip/N1.Lambda_tip_max ));                % Vesicles entering pool Neuron 1
    N2.Lambda_tip = N2.Lambda_tip ...                                           % Tip of Neuron 2
        + N1.dt*K*(- N2.alpha_r*(N2.Lambda_tip/N1.Lambda_tip_max)*(1-N2r1(end)-N2a1(end)) ...                 % Vesicles exiting pool Neuron 2
            + N2.beta_a*N2a1(end)*(1 - N2.Lambda_tip/N2.Lambda_tip_max));  

%     % Calculate total mass
%     Development_MassWholeSystem(i) = N1.Lambda_som + N1.Lambda_tip + N2.Lambda_tip  ... 
%         + (sum(N1r1(:)) + sum(N1a1(:)) + sum(N2r1(:)) + sum(N2a1(:)))*N2.dx; 
    
%     % Calculate entropy
%     eN1 = N1a1.*log(N1a1./N1.a_infinity) + N1r1.*log(N1r1./N1.r_infinity) + (1 - N1p1).*log((1 - N1p1)./(1-N1.a_infinity - N1.r_infinity)) ;
%     eN2 = N2a1.*log(N2a1./N2.a_infinity) + N2r1.*log(N2r1./N2.r_infinity) + (1 - N2p1).*log((1 - N2p1)./(1-N2.a_infinity - N2.r_infinity)) ;
%     Development_Entropy(i) = N1.dx*sum(eN1) + N2.dx*sum(eN2)  ...
%         + N1.Lambda_som*log(N1.Lambda_som/N1.Lambda_som_infinity) ...
%         + N1.Lambda_tip*log(N1.Lambda_tip/N1.Lambda_tip_infinity) ...
%         + N2.Lambda_tip*log(N2.Lambda_tip/N2.Lambda_tip_infinity);
   
end

toc

%% Test um Danilas Plots nachzubauen
ParticlesN1 = zeros(2,1);
ParticlesN2 = ParticlesN1;

ParticlesN1(1) = sum(N1a1(:))*N1.dx;
ParticlesN1(2) = sum(N1r1(:))*N1.dx;
ParticlesN2(1) = sum(N2a1(:))*N2.dx;
ParticlesN2(2) = sum(N2r1(:))*N2.dx;

figure(2);
plotN1a1 = subplot(1,4,1);
bar(ParticlesN1(1));
ylabel('Concentration of vesicles');
title('ANT N1')
ylim([0 1]);
plotN1a1.PlotBoxAspectRatio(2) = 2.5;
set(gca,'FontSize',15,'FontWeight','bold')

plotN2r1 = subplot(1,4,2);
bar(ParticlesN1(2));
title('RET N1')
ylim([0 1]);
plotN2r1.PlotBoxAspectRatio(2) = 2.5;
set(gca,'FontSize',15,'FontWeight','bold')

plotN1a1 = subplot(1,4,3);
bar(ParticlesN2(1));
title('ANT N2')
ylim([0 1]);
plotN1a1.PlotBoxAspectRatio(2) = 2.5;
set(gca,'FontSize',15,'FontWeight','bold')

plotN2r1 = subplot(1,4,4);
bar(ParticlesN2(2));
title('RET N2')
ylim([0 1]);
plotN2r1.PlotBoxAspectRatio(2) = 2.5;
set(gca,'FontSize',15,'FontWeight','bold')

% exportFigure('Danilasbildernachbauen2B.pdf', gcf);

%% Plot Pool Density Evolution in the pools
figure(3);
subplot(1,3,1);
plot(tt, Development_Lambda_som);
axis square
title('$\Lambda_{som}$');
xlabel('Time $t$');
ylabel('Concentration');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultAxesFontName', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(gca,'FontSize',12,'FontWeight','bold');

subplot(1,3,2);
plot(tt, Development_Lambda_tipN1);
axis square
title('$\Lambda_{N1}$');
xlabel('Time $t$');
ylabel('Concentration');
set(gca,'FontSize',12,'FontWeight','bold');

subplot(1,3,3),
plot(tt, Development_Lambda_tipN2);
axis square
title('$\Lambda_{N2}$');
xlabel('Time $t$');
ylabel('Concentration');
set(gca,'FontSize',12,'FontWeight','bold');

resizeFigure(gcf, [270,850]);
%exportFigure('Images/DevelopmentPoolConcentration2B.pdf', gcf);
disp('A figure of the Development of the Concentration in the pools was printed.')

%% Plot Development of Mass and Entropy of the system

figure(4);
plot(tt, Development_MassWholeSystem);
title('Mass Development of the whole System');
xlabel('Time $t$');
ylabel('Total Mass');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultAxesFontName', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(gca,'FontSize',20);
set(gca,'FontSize',17,'FontWeight','bold');

%saveas(gcf, 'Images/DevelopmentMass2B', 'pdf');
disp('A figure of the Development of the mass was printed.')

% figure(5);
% plot(tt, Development_Entropy);
% l = title('Relative Entropy');
% xl = xlabel('Time $t$');
% yl = ylabel('Total Entropy');
% set(groot, 'DefaultTextInterpreter', 'latex');
% set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
% set(groot, 'DefaultAxesFontName', 'latex');
% set(groot, 'DefaultLegendInterpreter', 'latex');
% set(l,'FontSize',20);
% set(xl,'FontSize',20);
% set(yl,'FontSize',20);
% set(gca,'FontSize',17,'FontWeight','bold');
% 
% %saveas(gcf, 'Images/DevelopmentEntropy2B', 'pdf');
% disp('A figure of the Development of the entropy was printed.')

