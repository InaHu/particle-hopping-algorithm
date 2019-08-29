function [a1, r1] = UpdateConcentrationParticleHopping(N) 
    % This function updates the concentration of anterograde and retrograde
    % moving particles with a partical hopping ansatz.
    %% Initialise
    N.p0 = N.a0 + N.r0;                         % Total concentration in the neurite
    
    % Initialise Additional Parameters for Scaling and Numerical Scheme
    C_a = 1/(2*N.eps_a);
    C_r = 1/(2*N.eps_r);
    tildeC = 1/N.typicalDiffusion;              
    tildeV = N.typicalPotential;
    
    % Initialse Vectors for new Concentrations
    a1 = zeros(length(N.a0),1); 
    r1 = zeros(length(N.r0),1);
    
    %% Solve in the Domain 
    % the parameter H and the function G are explained in the paper in Section 5
    H_a = 2*N.lambda_eps*(N.eps_a/(N.h.^2));         
    H_r = 2*N.lambda_eps*(N.eps_r/(N.h.^2));
    
    % Calculate new Concentrations
    G = - N.a0(2:end-1).*(1-N.p0(1:end-2)).*exp(- tildeC*tildeV*C_a*(N.V_a(2:end-1)-N.V_a(1:end-2)));          % Probability that an ANT-particle jumps away to the left
    G = G + N.a0(1:end-2).*(1-N.p0(2:end-1)).*exp(-tildeC*tildeV*C_a*(N.V_a(1:end-2)-N.V_a(2:end-1)));
    G = G - N.a0(2:end-1).*(1-N.p0(3:end)).*exp(-tildeC*tildeV*C_a*(N.V_a(2:end-1)-N.V_a(3:end)));
    G = G + N.a0(3:end).*(1-N.p0(2:end-1)).*exp(-tildeC*tildeV*C_a*(N.V_a(3:end)-N.V_a(2:end-1)));
    
    a1(2:end-1) = N.a0(2:end-1) + H_a*N.tau*G;
    
    G = - N.r0(2:end-1).*(1-N.p0(1:end-2)).*exp(-tildeC*tildeV*C_r*(N.V_r(2:end-1)-N.V_r(1:end-2)));          % Probability that a RET particle jumps away to the left
    G = G + N.r0(1:end-2).*(1-N.p0(2:end-1)).*exp(-tildeC*tildeV*C_r*(N.V_r(1:end-2)-N.V_r(2:end-1)));
    G = G - N.r0(2:end-1).*(1-N.p0(3:end)).*exp(-tildeC*tildeV*C_r*(N.V_r(2:end-1)-N.V_r(3:end)));
    G = G + N.r0(3:end).*(1-N.p0(2:end-1)).*exp(-tildeC*tildeV*C_r*(N.V_r(3:end)-N.V_r(2:end-1)));
    
    r1(2:end-1) = N.r0(2:end-1) + H_r*N.tau*G;

    %% Solve the Boundary Conditions

    % Boundary Flow at x = 0    
    a1(1) = N.a0(1) ...
        + N.tau*H_a*(-N.a0(1)*(1-N.p0(2))*exp(-tildeC*tildeV*C_a*(N.V_a(1)-N.V_a(2))) ...       
                    + N.a0(2)*(1-N.p0(1))*exp(-tildeC*tildeV*C_a*(N.V_a(2)-N.V_a(1)))) ... 
        + N.tau*N.lambda_in/N.h*N.alpha_a*((N.Lambda_som/N.Lambda_som_max)*(1-N.p0(1)));
    
    r1(1) = N.r0(1) ...
        + N.tau*H_r*(-N.r0(1)*(1-N.p0(2))*exp(-tildeC*tildeV*C_r*(N.V_r(1)-N.V_r(2))) ...
                    + N.r0(2)*(1-N.p0(1))*exp(-tildeC*tildeV*C_r*(N.V_r(2)-N.V_r(1)))) ...
        - N.tau*N.lambda_out/N.h*(N.beta_r*N.r0(1)*(1 - N.Lambda_som/N.Lambda_som_max));
    
    % Boundary Flow at x = 1 
    a1(end) = N.a0(end) ...
        + N.tau*H_a*(-N.a0(end)*(1-N.p0(end-1))*exp(-tildeC*tildeV*C_a*(N.V_a(end)-N.V_a(end-1))) ...
            + N.a0(end-1)*(1-N.p0(end))*exp(-tildeC*tildeV*C_a*(N.V_a(end-1)-N.V_a(end)))) ... 
        - N.tau*N.lambda_out/N.h*(N.beta_a*N.a0(end)*(1 - N.Lambda_tip/N.Lambda_tip_max)); 
    
    r1(end) = N.r0(end) ...
        + N.tau*H_r*(-N.r0(end)*(1-N.p0(end-1))*exp(-tildeC*tildeV*C_r*(N.V_r(end)- N.V_r(end-1))) ... 
            + N.r0(end-1)*(1-N.p0(end))*exp(-tildeC*tildeV*C_r*(N.V_r(end-1)-N.V_r(end)))) ... 
        + N.tau*N.lambda_in/N.h*(N.alpha_r*(N.Lambda_tip/N.Lambda_tip_max)*(1-N.p0(end)));
end