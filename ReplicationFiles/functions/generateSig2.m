function [sig2, returns, eps] = generateSig2(type, T, eps)

    if nargin < 3
        eps = randn(T,1);
    end
    
    switch type
       
        case 'constant'
            sig2 = 1;   
        
        case 'garch'
            omega = 0.05;
            alpha = 0.15;
            beta = 0.82;
            [returns, sig2] = generateGARCHVolatilities(T, omega, alpha, beta, eps);
    
        case 'cos'
            % entire wave length
            theta = 0.5;
            sig2 = 0.5+ theta*(1+cos(2*pi.*(1:T)'./T)).^2;

            % old parametrization sig2 = (1.5+cos((1:T)'/T*2*pi)).^2;
        
        case 'cos-down' % downward shicfting variance with smooth cosine function
            % half a wave length
            theta = 0.5;
            sig2 =  0.5 + theta*(1+cos(pi.*(1:T)'./T)).^2;
            sig2_print= 'cos(1.5+\frac{t}{2T\pi})';
     
        case 'cos-up' % upward shicfting variance with smooth cosine function
            % theta = 0.5;
            % sig2 = 0.5 +  theta*(1+ cos(pi+pi.*(1:T)'./T)).^2;
            theta = 0.5;
            sig2 = 2.5 - theta*(1+cos(2*pi.*(1:T)'./T)).^2;
    
        case 'ST-down'
            a       = 2;  % total change of the variance
            theta   = -0.25; % shift down or up
            lam     = 0.5;  % Date around which the change happens
            sig2 = 0.5+ a./(1+exp(-theta*((1:T)'-lam*T)));
        
        case 'ST-up'
            a       = 2;  % total change of the variance
            theta   = 0.25; % shift down or up
            lam     = 0.5;  % Date around which the change happens
            sig2 = 0.5+ a./(1+exp(-theta*((1:T)'-lam*T)));
    end

returns = sqrt(sig2).*eps;
end