function [J, Jq, Ts, Pv, B]=modelvmd0_janus(dm, por, r, Th, Pvc, h, alpha, ds)
    % sigma = alpha;
    % alpha = 0;
    sigma=0;

    J = zeros(length(dm),1);
    Jq = zeros(length(dm),1);
    Ts = zeros(length(dm),1);
    Pv = zeros(length(dm),1);
    B = zeros(length(dm),1);
    
    Th = Th+273; % K, bulk T
    
    hf = 1/(1/h+ds/0.2); % W/m2/K
    Hv = 2265*1e6/1e3/3600; % 629.1667 Wh/L, heat of evaporation
    K = 1/(1/150+ds/1e-9/3600/1e3); % LMH, mass transfer coeff
    Cb = 0.6; % molality of hot side: mol NaCl per kg water

    Pvs = @(x) exp(23.5377-4016.3622/(373-38.6339))*exp(40650/8.314/373)*exp(-40650*(1-alpha)/8.314./x);  % Pa, saturated vapor pressure of water
    %Pvs = @(x) exp(23.5377-4016.3622/(373-38.6339))*exp(-40650/8.314*(1./x-1/373));  % Pa, saturated vapor pressure of water
    
    aw = @(x) 1-0.03112*x-0.001482*(x)^2;
    eta = @(x) 1 + (x^2/4) - (x/4)*(x^2+4)^0.5 - ((8-x^2)*(x^2+4)^0.5+x^3-16)^2/(72*x*(x^2+4)^0.5-288*log(x+(x^2+4)^0.5)+288*log(2));
    tor = por^-0.5;



    for i = 1:length(dm)

        % x1: Jw, LMH
        % x2: Jq, W/m2
        % x3: Tsh, K
        % x4: cm, mol/kg

        F = @(x) [

            x(2) - hf*(Th-abs(x(3)));

            x(2) - x(1)*Hv*(1-alpha);

            x(1) - (Pvs(abs(x(3)))*aw(abs(x(4)))-Pvc)*por*3600*0.018/((2*3.14*0.018*8.314*abs(x(3)))^0.5*(1/condcoeff(abs(x(3)))-1)*(1-sigma) +(2*3.14*0.018*8.314*abs(x(3)))^0.5*(1/eta(dm(i)*tor/r))  );

            abs(x(4)) - Cb*exp(x(1)/K);
        ];
        
        options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
        options.MaxIterations = 30000;
        options.MaxFunctionEvaluations = 1000000;
        options.StepTolerance = 1e-16;
        options.FunctionTolerance = 1e-16;
        
        if i == 1
            x0 = [0,0,Th-5,Cb];
        else
            x0 = x;
        end
        x = fsolve(F,x0,options); % solve systems of equations
        
        [i,sum(F(x).^2),condcoeff(abs(x(3)))]

        J(i) = x(1);
        Jq(i) = x(2);
        Ts(i) = abs(x(3))-273; % 
        Pv(i) = Pvs(abs(x(3)))*aw(abs(x(4))); % Pa

        P1 = Pvs(abs(x(3)))*aw(abs(x(4))) - x(1)*(2*3.14*0.018*8.314*abs(x(3)))^0.5*(1/condcoeff(abs(x(3)))-1)/por/3600/0.018;
        
        
        B(i) = x(1)/(Pvs(abs(x(3)))*aw(abs(x(4)))-Pvc); % LMH/Pa

    end

end