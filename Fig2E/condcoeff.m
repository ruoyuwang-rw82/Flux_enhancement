function y=condcoeff(T)
    T = T-273;
    if T<40
        y = -0.014132684*T+0.861518592;
    elseif T>50
        y = -0.00279134*T+0.356861903;
    else
        y = -0.007891635*T+0.611876656;
    end
    
   
    %y=0.5;
end