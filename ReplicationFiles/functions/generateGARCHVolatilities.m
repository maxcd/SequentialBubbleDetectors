function [returns, vol] = generateGARCHVolatilities(T, omega, alpha, beta, eps)
    
    % Function to generate GARCH(1,1) volatilities
    % T: Number of periods to simulate
    % omega:Constant term in GARCH model
    %       Generally, values range from 0.1 to 0.3. This parameter captures the 
    %       impact of previous shocks (squared returns) on current volatility.
   
    % alpha:Coefficient for lagged squared returns
    %       Generally, values range from 0.1 to 0.3. This parameter captures the 
    %       impact of previous shocks (squared returns) on current volatility.
    % beta: Coefficient for lagged conditional variance
    %       Commonly between 0.6 and 0.9. This reflects the influence of past 
    %       conditional variances on current volatility.

    % Initialize variables
    returns = zeros(T, 1);
    h = zeros(T, 1); % Conditional variance

    % Initial conditions
    h(1) = omega / (1 - alpha - beta); % Starting value for the variance
    returns(1) = sqrt(h(1)) .* eps(1); 

    for t = 2:T
        % Update the conditional variance
        h(t) = omega + alpha * returns(t-1).^2 + beta * h(t-1);
        % Generate returns based on volatility
        returns(t) = sqrt(h(t)) .* eps(t); 
    end
    vol = h;
end