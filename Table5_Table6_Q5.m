function [] = Table5_Table6_Q5( mydigits )
 
clc;
fprintf( 'You have set the precision as %d digits.\r', mydigits );
 
close all;
OutPutPath = '.\';

warning( 'on', 'MATLAB:nearlySingularMatrix' );
warning( 'off', 'MATLAB:mpath:nameNonexistentOrNotADirectory' );

ToolboxPath = 'D:\Multiprecision Computing Toolbox\';
addpath( ToolboxPath );
[msgstr, msgid] = lastwarn;
if strcmp(msgid, 'MATLAB:mpath:nameNonexistentOrNotADirectory')
    fprintf('ToolboxPath Wrong!!!\rPleae specify the right path that installs the "Multiprecision Computing Toolbox".\r');
    return;
end 

for kkk = 1
    try
        mp.Digits( mydigits );
    catch ME
        if ~isempty( strfind(ME.message, 'expired') )
            disp(ME.message);
        else 
            fprintf( ['You need first install the "Multiprecision Computing Toolbox."\r', ...
                'Please you can download the toolbox from "https://www.advanpix.com/"\r'] );
        end
        return;
    end
    
    ImUnit = mp( sqrt(-1) );
    x0 = mp( 0 );
    mu = mp( 0 );
    sigma = 3/mp( 10 );
    lambda = mp( 5 );
    hitlevel = 1/mp( 10 ); % hitting level
    
    eta = mp( linspace( mp(20), mp(40), 10 )' ); % exp parameters for upward jumps 
    p = mp( ones(size(eta)) ) / length( eta ) / 2; % probability for upward jumps
    theta = eta; % exp parameters for downward jumps
    q = p; % probability for downward jumps
     
    
    mubar = mu + lambda * (p'*(1./eta) - q'*(1./theta));
    zeta = sigma^2/2 + lambda * (p'*(1./eta.^2) + q'*(1./theta.^2));
    varpi = lambda * (p'*(1./eta.^3) - q'*(1./theta.^3)); 
    
    length_eta = mp( length( eta ) );
    length_theta = mp( length( theta ) );
    length_coeff = mp( length_eta + length_theta + 2 + 1 );
    nroots = mp( length_eta + length_theta + 2 );
    
    y = 1/mp( 10 ); % overshot level;
    alpha = mp( min(eta) ) / 2; % For J_4 ( Laplace Transform Parameter )
    
    br = [ [       q*0;        ones(size(p));   0;   1 ], ...
        [          q*0;         exp(-eta.*y);   0;   0 ], ...
        [    -1./theta;                  p*0;  -1;   0 ], ...
        [          q*0;   eta./(eta - alpha);   0;   1 ], ...
        [          q*0;               1./eta;   0;   0 ] ];
    br = mp( br );
    d0 = mp( [0./theta; eta./eta; 0; 1]);
    d1 = mp( [1./theta; 1./eta; 1; 0]);
    d2 = mp( [-2*hitlevel./theta - 2./theta.^2; 2./eta.^2; -2*hitlevel; 0]);
    d3 = mp( [3*hitlevel.^2./theta + 6*hitlevel./theta.^2 + 6./theta.^3; 6./eta.^3; 3*hitlevel.^2; 0]);
    
    onesn1 = mp( ones(nroots, 1) );
    onesn2 = mp( ones(nroots, 1) );
    onesn1( length_theta + 1 ) = 0;
    onesn2( length_theta + 2 ) = 0;
    onesn1 = mp ( onesn1 );
    onesn2 = mp ( onesn2 );
    
    mus = mp( [0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 1e-5, 1e-6, 1e-7, 1e-8, 0] );
    betas = mp( [0.01, 0.001, 0.0005, 0.0002, 0.0001, 1e-5, 1e-6, 1e-7, 1e-8, 0] );
    mus = 1./mp(round(1./mus));
    betas = 1./mp(round(1./betas));
     
    length_mu = mp( length(mus) );
    length_beta = mp( length(betas) );
    
    J = mp( zeros(1, 5) );
    Etau = mp( zeros(length_mu, length_beta) );
    Qv = mp( zeros(length_mu, length_beta) );
    
    tic
    
    for mu_th = 1:length(mus)
        for beta_th = 1:length(betas)
            mu = mp( mus(mu_th) );
            beta = mp( betas(beta_th) );
            mubar = mp( mu + lambda * (p'*(1./eta) - q'*(1./theta)) );
            
            LundbergCoeff = PowerExpandForLundbergEq( mu, sigma, lambda, eta, p, theta, q, beta, length_eta, length_theta, length_coeff );
            Lundbergroots = roots( LundbergCoeff );
            [~, indextemp] = sort( real(Lundbergroots) );
            Lundbergroots = Lundbergroots( indextemp );
            
            if beta==0
                if mubar>0
                    Lundbergroots(length_theta + 2) = 0;
                elseif mubar<0
                    Lundbergroots(length_theta + 1) = 0;
                elseif mubar==0
                    Lundbergroots(length_theta + 1:length_theta + 2) = 0;
                end
            end
            
            for i = nroots:-1:1
                A(:, i) = mp( h( Lundbergroots(i), theta, eta, hitlevel ) );
            end
            
            if beta~=0
                C = A\br;
                J = exp( - ( hitlevel - x0 ) * Lundbergroots' ) * C; 
                Qv( mu_th, beta_th ) = ( beta * ( x0 - hitlevel*J(1) + J(3) - J(5) ) + mubar * ( mp(1) - J(1) ) )./beta.^2;
            else
                [ D00, D01, D02, D03, D10, D20, D12, ...
                    D10ik, D20ik, D01ik, D02ik, D03ik, D12ik, D20d1k ] = DA( A, d1, d2, d3, br, nroots, length_theta );
                
                ExpRoots_x0 = exp( - ( hitlevel - x0 ) * Lundbergroots' );
                ExpRoots_x0_n1 = ExpRoots_x0 .* onesn1';
                ExpRoots_x0_n2 = ExpRoots_x0 .* onesn2';
                
                if mubar==0
                    J = - ( hitlevel - x0 ) * D10ik( length_theta + 1, : )/D10 ...
                        + ExpRoots_x0_n1 * D10ik/D10;
                    Etau( mu_th, beta_th ) = -1/2./zeta .* ( ( hitlevel - x0 ).^2 + ( hitlevel - x0 ) * D20/D10 ...
                        + ExpRoots_x0_n1 * D20d1k/D10 ); 
                elseif mubar~=0
                    for i = nroots:-1:1
                        A1(:, i) = h1( Lundbergroots(i), mu, sigma, lambda, eta, p, theta, q, hitlevel );
                        delta1(i,:) = 1/G1( Lundbergroots(i), mu, sigma, lambda, eta, p, theta, q );
                    end
                    C = A\br;
                    C1 = A\( - A1 * C(:,1) );
                    J = ExpRoots_x0 * C;
                    Etau( mu_th, beta_th ) = ExpRoots_x0 * ( ( hitlevel - x0 ) .* delta1 .* C(:,1) - C1 ); 
                end
                
                if mubar>0
                    Qv( mu_th, beta_th ) = 1./mubar/2.*( hitlevel.^2 - x0.^2 + ( 2*hitlevel*D01 + D02 )/D00 ...
                        - ExpRoots_x0_n2 * ( 2*hitlevel.*D01ik(:,1) + D02ik(:,1) )/D00 ) ...
                        + zeta.*mp( ( x0 - hitlevel + J(3) - J(5) )./mubar.^2 );
                elseif mubar==0
                    Qv( mu_th, beta_th ) = 1./zeta/2.*( varpi./zeta.*( ( hitlevel - x0 ).^2 + ( hitlevel - x0 ).*D02/D01 ...
                        - ExpRoots_x0_n2 * ( D02ik(:,3) - D02ik(:,5) )/D01 ) ...
                        - ( hitlevel - x0 ).*D03./D01/3 - hitlevel.*( hitlevel - x0 + mp(1) ).*D02/D01 ...
                        - ( hitlevel - x0 ).^2.*( x0 + 2*hitlevel )/3 - hitlevel.*D12/D01 ...
                        + ExpRoots_x0_n2 * ( 3.*hitlevel.*D12ik(:,1) + D03ik(:,3) - D03ik(:,5) )./D01/3 );
                elseif mubar<0
                    Qv( mu_th, beta_th ) = 1./mubar/2.*( hitlevel.^2 - x0.^2 + ( 2*hitlevel*D10 + D20 )/D00 ...
                        - ExpRoots_x0_n1 * ( 2*hitlevel.*D10ik(:,1) + D20ik(:,1) )/D00 ) ...
                        + zeta./mubar.^2.*( x0 - hitlevel + J(3) - J(5) );
                end  
            end 
        end
    end
    
    CPUtime = toc; 
    TureValue = round( Qv*1e6 )/1e6;  
    TureValue = double(TureValue); 
end
fprintf( 'Values with %d digits precision are as follows:\r\r', mydigits );
disp( num2str( TureValue, '%13.4e' ) );
return;
 
end
  
%%%%%%%%%%%%%%%%%%%% Function h ( Columns of Matrix A(beta) ) %%%%%%%%%%%%%%%%%%%%%%%%
function [ hx ] = h( x, theta, eta, hitlevel )
hx = mp( [ [ x ./ ( theta + x ) .* exp( - x .* hitlevel ) ]; ...
    [                              eta ./ ( eta - x ) ]; ...
    [                     x .* exp( - x .* hitlevel ) ]; ...
    [                                               1 ];     ] ); %#ok<*NBRAK>
end

%%%%%%%%%%%%%%%%%%%% Function h1 ( The first order derivative of the Columns of Matrix A(beta) ) %%%%%%%%%%%%%%%%%%%%%%%%
function [ hx ] = h1( x, mu, sigma, lambda, eta, p, theta, q, hitlevel )
hx = [   [              - x .* exp( - x .* hitlevel ) ./ ( theta + x ).^2  ...
    +                          exp( - x .* hitlevel ) ./ ( theta + x ) ...
    -         hitlevel .* x .* exp( - x .* hitlevel ) ./ ( theta + x )  ]; ...
    [                                             eta ./ ( eta - x ).^2 ]; ...
    [  exp( - x .* hitlevel ) - hitlevel .* x .* exp( - x .* hitlevel ) ]; ...
    [                                                                 0 ];     ] ...
    /G1( x, mu, sigma, lambda, eta, p, theta, q ); %#ok<*NBRAK>
end

%%%%%%%%%%%%%%%%%%%% Function G ( Lundberg equation ) %%%%%%%%%%%%%%%%%%%%%%%%
function [ Gx ] = G( x, mu, sigma, lambda, eta, p, theta, q )
Gx = mu.*x + 1/2*sigma.^2.*x.^2 + lambda.*( p' * ( eta./(eta-x) ) + q' * ( theta./(theta+x) ) - 1 );
end

%%%%%%%%%%%%%%%%%%%% Function G1 ( the first order derivative of the Lundberg equation ) %%%%%%%%%%%%%%%%%%%%%%%%
function [ Gx ] = G1( x, mu, sigma, lambda, eta, p, theta, q )
Gx = mu + sigma.^2.*x + lambda.*( p' * ( eta./(eta-x).^2 ) - q' * ( theta./(theta+x).^2 ) );
end

%%%%%%%%%%%%%%%%%%%% Find determinants of various matrix %%%%%%%%%%%%%%%%
function [ D00, D01, D02, D03, D10, D20, D12, ...
    D10ik, D20ik, D01ik, D02ik, D03ik, D12ik, D20d1k ] = DA( A, d1, d2, d3, br, nroots, length_theta )
A01 = A;
A02 = A;
A03 = A;
A10 = A;
A20 = A;
A12 = A;

A01( :, length_theta + 2 ) = d1;
A02( :, length_theta + 2 ) = d2;
A03( :, length_theta + 2 ) = d3;
A10( :, length_theta + 1 ) = d1;
A20( :, length_theta + 1 ) = d2;
A12( :, length_theta + 1 ) = d1;
A12( :, length_theta + 2 ) = d2;

D00 = det(A);

D01 = det(A01);
D02 = det(A02);
D03 = det(A03);
D10 = det(A10);
D20 = det(A20);
D12 = det(A12);

D10ik = zeros(nroots, 5);
D20ik = zeros(nroots, 5);
D01ik = zeros(nroots, 5);
D02ik = zeros(nroots, 5);
D03ik = zeros(nroots, 5);
D12ik = zeros(nroots, 5);
D20d1k = zeros(nroots, 1);

for k=1:nroots
    for i=1:5
        A10ik = A10;
        A20ik = A20;
        A01ik = A01;
        A02ik = A02;
        A03ik = A03;
        A12ik = A12;
        
        A10ik(:,k) = br(:,i);
        A20ik(:,k) = br(:,i);
        A01ik(:,k) = br(:,i);
        A02ik(:,k) = br(:,i);
        A03ik(:,k) = br(:,i);
        A12ik(:,k) = br(:,i);
        
        D10ik(k, i) = det(A10ik);
        D20ik(k, i) = det(A20ik);
        D01ik(k, i) = det(A01ik);
        D02ik(k, i) = det(A02ik);
        D03ik(k, i) = det(A03ik);
        D12ik(k, i) = det(A12ik);
    end
    A20d1k = A20;
    A20d1k(:,k) = d1;
    D20d1k(k,:) = det(A20d1k);
end
end


%%%%%%%%%%%%%%%%%%%% Find coeffecients of Lundberg equation %%%%%%%%%%%%%%%%
function [ LundbergCoeff ] = PowerExpandForLundbergEq( mu, sigma, lambda, eta, p, theta, q, beta, length_eta, length_theta, length_coeff )
%  %%%%%%%%%%%%%% Input %%%%%%%%%%%%%%
%    eta: exp parameters for upward jumps
%      p: probability for upward jumps
%  theta: exp parameters for downward jumps
%      q: probability for downward jumps 
%     mu: drift parameter
%  sigma: volatility parameter
% lambda: jump intensity
%   beta: discount factor  
%  %%%%%%%%%%%%%% Output %%%%%%%%%%%%%%
%  LundbergCoeff: Coefficients of Lundberg Equation
%  containing length(eta)+length(theta)+2+1 polynomial coefficients, 
%  starting with the coefficient of x^(length(eta)+length(theta)+2). 

Coeff = zeros( length_coeff, 3 );
Coeff(1:3, 1) = [ - beta - lambda; mu; sigma^2/2 ];
Coeff(1, 2) = 1;
Coeff(1:3, 3) = [ - beta - lambda; mu; sigma^2/2 ];
for i = 1:length_eta
    Coeff( :, 3 ) = ( - 1 ) * [0; Coeff( 1:end-1, 1 )] + eta(i) * Coeff( :, 1 ) + lambda * p(i) * eta(i) * Coeff(:, 2);
    Coeff( :, 1 ) = Coeff( :, 3 );
    Coeff( :, 2 ) = ( - 1 ) * [0; Coeff( 1:end-1, 2 )] + eta(i) * Coeff( :, 2 );
end

for i = 1:length_theta
    Coeff( :, 3 ) = [0; Coeff( 1:end-1, 1 )] + theta(i) * Coeff( :, 1 ) + lambda * q(i) * theta(i) * Coeff(:, 2);
    Coeff( :, 1 ) = Coeff( :, 3 );
    Coeff( :, 2 ) = [0; Coeff( 1:end-1, 2 )] + theta(i) * Coeff( :, 2 );
end

LundbergCoeff = flipud( Coeff( :, 3 ) );

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
