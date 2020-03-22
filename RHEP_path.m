function [] = RHEP_path(  )

clc; 
close all;

x0 = 0;
mu = 0;
sigma = 0.16;
lambda = 5;  

eta = linspace( 20, 40, 10 )'; % exp parameters for upward jumps
p = ones(size(eta))/length( eta )/2; % probability for upward jumps
theta = eta; % exp parameters for downward jumps
q = p; % probability for downward jumps 

quantile_hep = Find_quantile_hep( eta, p, theta, q );

t_all = 1;
stepsize = 1/300000; 
ts = 0:stepsize:t_all;
Nstep = round( t_all/stepsize );
Nsample = 3;
Nexperiment = 1;

sigmadt = sigma * sqrt(stepsize); 

ujump = rand( Nstep, Nsample );
ujump = min( 1e6, floor( ujump*1e6 ) + 1 );
ujump = quantile_hep( ujump ); 

for i = 1:Nexperiment
    RandNormal = randn( Nstep, Nsample );
    RandUnit = rand( Nstep, Nsample );
    ulambda = RandUnit <= lambda * stepsize; 
    
    wn = sigmadt * RandNormal + mu * stepsize;
    
    Xt = cumsum( wn + ulambda.*ujump ); 
    Xt = x0 + [Xt(1,:)*0; Xt];  
    Lt = - min( 0, cummin(Xt) );
    hXt = Xt + Lt;
    plot(ts, hXt)
end 

return;

end

%%%%%% Find Quantiles for HEPs' Jump Distribution ( To Speed up the Simulation )
function [ quantile_hep ] = Find_quantile_hep( eta, p, theta, q )
% cdf_hyper_exp( -0.656236, eta, p, theta, q );
myquantile = [ 1e-10; 2e-10; 5e-10; 1e-9; 2e-9; 5e-9; 1e-8; 2e-8; 5e-8; ...
    1e-7; 2e-7; 5e-7; 1e-6; 2e-6; 5e-6; (1e-5:1e-5:9e-5)'; (1e-4:1e-4:9e-4)'; ...
    (1e-3:1e-3:1-1e-3)';
    (1-9e-4:1e-4:1-1e-4)'; sort(1- [1e-10; 2e-10; 5e-10; 1e-9; 2e-9; 5e-9; 1e-8; 2e-8; 5e-8; ...
    1e-7; 2e-7; 5e-7; 1e-6; 2e-6; 5e-6; (1e-5:1e-5:9e-5)']); ] ;
%     myquantile_full = [(1e-10:1e-10:9e-10)'; (1e-9:1e-9:9e-9)'; (1e-8:1e-8:9e-8)'; ...
%         (1e-7:1e-7:9e-7)'; (1e-6:1e-6:1-1e-6)'; (1-9e-7:1e-7:1-1e-7)'; ...
%         (1-9e-8:1e-8:1-1e-8)'; (1-9e-9:1e-9:1-1e-9)'; (1-9e-10:1e-10:1-1e-10)'; ];
myquantile_full = [1e-7; (1e-6:1e-6:1-1e-6)'; ];
quantile_hep = myquantile*0;
x0 = -0.656236;
for i=1:length(myquantile)
    quantile_hep(i) = fzero( @(x)( cdf_hyper_exp( x, eta, p, theta, q ) - myquantile(i) ), x0 );
    x0 = quantile_hep(i);
    if mod(i, 1000) == 0
        fprintf( '[%0.1e, %0.6f]\r', [i, quantile_hep(i)] );
    end
end
quantile_hep = pchip( myquantile, quantile_hep, myquantile_full);

end

%%%%%%%%%%%%%%%%%%%% Function h ( Columns of Matrix A(beta) ) %%%%%%%%%%%%%%%%%%%%%%%%
function [ cdf ] = cdf_hyper_exp( x, eta, p, theta, q )

if x<=0
    cdf = q' * exp( theta * x );
else
    cdf = sum(q) + p' * ( 1 - exp( - eta * x ) );
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
