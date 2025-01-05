function [J, Jq, Ts, Pv, B]=model0_janus(dm, por, r, kp, Th, Tc, h, alpha, ds, face)
    % sigma = alpha;
    % alpha = 0;
    sigma=0;

    J = zeros(length(dm),1);
    Jq = zeros(length(dm),1);
    Ts = zeros(length(dm),2);
    Pv = zeros(length(dm),4);
    B = zeros(length(dm),1);
    
    Th = Th+273; % K, bulk T
    Tc = Tc+273;
    
    if face == 'F'
        hf = 1/(1/h+ds/0.2); % W/m2/K
        hp = h; 
        K = 1/(1/150+ds/1e-9/3600/1e3); % LMH, mass transfer coeff
        betaf = 0;
        betad = alpha;
        Pvs = @(x) exp(23.5377-4016.3622/(373-38.6339))*exp(-40650/8.314*(1./x-1/373));  % Pa, saturated vapor pressure of water
        Pvsh = @(x) exp(23.5377-4016.3622/(373-38.6339))*exp(40650/8.314/373)*exp(-40650*(1-alpha)/8.314./x);  % Pa, saturated vapor pressure of water
        %Pvsh = Pvs;

    else
        hf = h; % W/m2/K
        hp = 1/(1/h+ds/0.2); 
        K = 150; % LMH, mass transfer coeff
        betaf = alpha;
        betad = 0;
        Pvsh = @(x) exp(23.5377-4016.3622/(373-38.6339))*exp(-40650/8.314*(1./x-1/373));  % Pa, saturated vapor pressure of water
        Pvs = @(x) exp(23.5377-4016.3622/(373-38.6339))*exp(40650/8.314/373)*exp(-40650*(1-alpha)/8.314./x);  % Pa, saturated vapor pressure of water
        %Pvsh = Pvs;

    end

    km = por*0.02735+(1-por)*kp; % W/m/K, thermal conductivity
    hm = km./dm; % W/m2/K
    Hv = 2256*1e6/1e3/3600; % 629.1667 Wh/L, heat of evaporation
    
    Cb = 0.6; % molality of hot side: mol NaCl per kg water
    

    aw = @(x) 1-0.03112*x-0.001482*(x)^2;
    eta = @(x) 1 + (x^2/4) - (x/4)*(x^2+4)^0.5 - ((8-x^2)*(x^2+4)^0.5+x^3-16)^2/(72*x*(x^2+4)^0.5-288*log(x+(x^2+4)^0.5)+288*log(2));
    tor = por^-0.5;



    for i = 1:length(dm)

        % x1: Jw, LMH
        % x2: Jq, W/m2
        % x3: Tsh, K
        % x4: Tsc, K
        % x5: cm, mol/kg

        F = @(x) [

            x(2) - hf*(Th-abs(x(3))) + x(1)*Hv*betaf;

            x(2) - hp*(abs(x(4))-Tc) + x(1)*Hv*betad;

            x(2) - x(1)*Hv*(1-alpha) - hm(i)*(abs(x(3))-abs(x(4)));

            x(1) - (Pvsh(abs(x(3)))*aw(abs(x(5)))-Pvs(abs(x(4))))*por*3600*0.018/((2*3.14*0.018*8.314*(abs(x(3))+abs(x(4)))/2)^0.5*(1/condcoeff(abs(x(3)))+1/condcoeff(abs(x(4)))-2)*(1-sigma) + (8.314*(abs(x(3))+abs(x(4)))/2*dm(i)*tor*(1-(Pvsh(abs(x(3)))*aw(abs(x(5)))+Pvs(abs(x(4))))/2/1e5))/(1.87e-10*((abs(x(3))+abs(x(4)))/2)^2.07) +(2*3.14*0.018*8.314*(abs(x(3))+abs(x(4)))/2)^0.5*(1/eta(dm(i)*tor/r))  );

            abs(x(5)) - Cb*exp(x(1)/K);
        ];
        
        options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
        options.MaxIterations = 30000;
        options.MaxFunctionEvaluations = 1000000;
        options.StepTolerance = 1e-16;
        options.FunctionTolerance = 1e-16;
        
        if i == 1
            x0 = [0,0,Th-5,Tc+5,Cb];
        else
            x0 = x;
        end
        x = fsolve(F,x0,options); % solve systems of equations
        
        [i,sum(F(x).^2),condcoeff(abs(x(3))),condcoeff(abs(x(4)))]

        J(i) = x(1);
        Jq(i) = x(2);
        Ts(i,1:2) = [abs(x(3)), abs(x(4))]-273; % 

        P1 = Pvsh(abs(x(3)))*aw(abs(x(5))) - x(1)*(2*3.14*0.018*8.314*(abs(x(3))+abs(x(4)))/2)^0.5*(1/condcoeff(abs(x(3)))-1)/por/3600/0.018;
        P2 = Pvs(abs(x(4))) + x(1)*(2*3.14*0.018*8.314*(abs(x(3))+abs(x(4)))/2)^0.5*(1/condcoeff(abs(x(4)))-1)/por/3600/0.018;
        Pv(i,1:4) = [Pvsh(abs(x(3)))*aw(abs(x(5))),P1,P2, Pvs(abs(x(4)))]; % Pa
        
        B(i) = x(1)/(Pvsh(abs(x(3)))*aw(abs(x(5)))-Pvs(abs(x(4)))); % LMH/Pa

    end

end