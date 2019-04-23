function [a1, r1, j] = UpdateConcentrationParticleHopping(N) 
    % This function updates the concentration of anterograde and retrograde
    % moving particles with a partical hopping ansatz.
    %% Initialise
    N.p0 = N.a0 + N.r0;                         % Total concentration in the neurite
    
    % Initialise additional Parameters for Scaling and Numerical Scheme
    C_a = 1/(2*N.eps_a);
    C_r = 1/(2*N.eps_r);
    tildeC = 1/N.typicalDiffusion;
    
    % Initialse Vectors for new concentrations
    a1 = zeros(length(N.a0),1); 
    r1 = zeros(length(N.r0),1);
    
    %% Solve in the Domain
    H_a = (N.typicalTime*N.typicalDiffusion)/(2*N.typicalLength^2*C_a*N.dx^2);
    H_r = (N.typicalTime*N.typicalDiffusion)/(2*N.typicalLength^2*C_a*N.dx^2);
    
    % Calculate new concentrations
    G = - N.a0(2:end-1).*(1-N.p0(1:end-2)).*exp(- tildeC*N.typicalPotential*C_a*(N.V_a(2:end-1)-N.V_a(1:end-2)));          % Probability that an ANT-particle jumps away to the left
    G = G + N.a0(1:end-2).*(1-N.p0(2:end-1)).*exp(-tildeC*N.typicalPotential*C_a*(N.V_a(1:end-2)-N.V_a(2:end-1)));
    G = G - N.a0(2:end-1).*(1-N.p0(3:end)).*exp(-tildeC*N.typicalPotential*C_a*(N.V_a(2:end-1)-N.V_a(3:end)));
    G = G + N.a0(3:end).*(1-N.p0(2:end-1)).*exp(-tildeC*N.typicalPotential*C_a*(N.V_a(3:end)-N.V_a(2:end-1)));
    
    a1(2:end-1) = N.a0(2:end-1) + H_a*N.dt*G;
    
    G = - N.r0(2:end-1).*(1-N.p0(1:end-2)).*exp(-tildeC*N.typicalPotential*C_r*(N.V_r(2:end-1)-N.V_r(1:end-2)));          % Probability that a RET particle jumps away to the left
    G = G + N.r0(1:end-2).*(1-N.p0(2:end-1)).*exp(-tildeC*N.typicalPotential*C_r*(N.V_r(1:end-2)-N.V_r(2:end-1)));
    G = G - N.r0(2:end-1).*(1-N.p0(3:end)).*exp(-tildeC*N.typicalPotential*C_r*(N.V_r(2:end-1)-N.V_r(3:end)));
    G = G + N.r0(3:end).*(1-N.p0(2:end-1)).*exp(-tildeC*N.typicalPotential*C_r*(N.V_r(3:end)-N.V_r(2:end-1)));
    
    r1(2:end-1) = N.r0(2:end-1) + H_r*N.dt*G;

    %% Solve the Boundary Conditions
    H_a = (N.typicalDiffusion)/(2*C_a*N.dx^2);
    H_r = (N.typicalDiffusion)/(2*C_a*N.dx^2);
    
    %lambda_B = (N.typicalTime*N.typicalInOutFlux)/(N.typicalLength*N.typicalConcentration);
    lambda_B = N.typicalInOutFlux/N.typicalConcentration;
  
    % Boundary Flow at x = 0    
    a1(1) = N.a0(1) ...
        + H_a*(-N.a0(1)*(1-N.p0(2))*exp(-tildeC*N.typicalPotential*C_a*(N.V_a(1)-N.V_a(2))) ... 
        	+ N.a0(2)*(1-N.p0(1))*exp(-tildeC*N.typicalPotential*C_a*(N.V_a(2)-N.V_a(1))) ... 
        	+ lambda_B*tildeC*C_a*N.dx*N.Lambda_som*N.alpha_a*(1-N.p0(1)));
    
    r1(1) = N.r0(1) ...
        + H_r*(-N.r0(1)*(1-N.p0(2))*exp(-tildeC*N.typicalPotential*C_r*(N.V_r(1)-N.V_r(2))) ...
            + N.r0(2)*(1-N.p0(1))*exp(-tildeC*N.typicalPotential*C_r*(N.V_r(2)-N.V_r(1))) ...
            - lambda_B*tildeC*C_r*N.beta_r*N.r0(1)*(N.Lambda_som_max - N.Lambda_som));
    
    % Boundary Flow at x = 1 
    a1(end) = N.a0(end) ...
        + H_a*(-N.a0(end)*(1-N.p0(end-1))*exp(-tildeC*N.typicalPotential*C_a*(N.V_a(end)-N.V_a(end-1))) ...
            + N.a0(end-1)*(1-N.p0(end))*exp(-tildeC*N.typicalPotential*C_a*(N.V_a(end-1)-N.V_a(end))) ... 
            - lambda_B*tildeC*C_a*N.beta_a*N.a0(end)*(N.Lambda_tip_max - N.Lambda_tip));
    
    r1(end) = N.r0(end) ...
        + H_r*(-N.r0(end)*(1-N.p0(end-1))*exp(-tildeC*N.typicalPotential*C_r*(N.V_r(end)- N.V_r(end-1))) ... 
            + N.r0(end-1)*(1-N.p0(end))*exp(-tildeC*N.typicalPotential*C_r*(N.V_r(end-1)-N.V_r(end))) ... 
            + lambda_B*tildeC*C_r*N.alpha_r*N.Lambda_tip*(1-N.p0(end)));
    
    j = N.a0 + N.r0 - a1 - r1;
    %disp(['What enters the neurite = ' num2str( N.eps_a*N.dt/(N.dx)^2*C_a*N.dx*N.Lambda_som*N.alpha_a*(1-N.p0(1))) ', (1-\rho)= ' num2str((1-N.p0(1)))]);
end