function [] = Figure1_Q2_to_Q5()
clc;
%% get the MCVs
MCtrue = dir;
MCtrue = ismember( 'QvEtau_MC.mat', {MCtrue.name}' );
if ~MCtrue
    RHEP_MC( );
end

close all;
chvar_name = 'beta'; % 'mu', 'beta'

OutPutPathFigure = '.\';
warning( 'off', 'MATLAB:nearlySingularMatrix' );
warning( 'off', 'MATLAB:mpath:nameNonexistentOrNotADirectory' );

%% get the LIVs
ToolboxPath = 'D:\Multiprecision Computing Toolbox\';
addpath( ToolboxPath );
[msgstr, msgid] = lastwarn;
if strcmp(msgid, 'MATLAB:mpath:nameNonexistentOrNotADirectory')
    fprintf('ToolboxPath Wrong!!!\rPleae specify the right path that installs the "Multiprecision Computing Toolbox".\r');
    return;
end

mydigits = 200;

tic;

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
x0 = mp(0);
mu = mp(0);
sigma = 3/mp(10);
lambda = mp(5);
hitlevel = 1/mp(10); % hitting level
y = 1/mp(100); % overshot level;

eta = mp( linspace( mp(20), mp(40), 10 )' ); % exp parameters for upward jumps
p = mp( ones(size(eta)) ) / length( eta ) / 2; % probability for upward jumps
theta = eta; % exp parameters for downward jumps
q = p; % probability for downward jumps
alpha = min(eta)/2; % For J_4 ( Laplace Transform Parameter )

mubar = mu + lambda * (p'*(1./eta) - q'*(1./theta));
zeta = sigma^2/2 + lambda * (p'*(1./eta.^2) + q'*(1./theta.^2));
varpi = lambda * (p'*(1./eta.^3) - q'*(1./theta.^3));
fprintf( 'mubar = %0.4f, zeta = %0.4f, varpi = %0.4f\r', [mubar, zeta, varpi] );

length_eta = mp( length( eta ) );
length_theta = mp( length( theta ) );
length_coeff = mp( length_eta + length_theta + 2 + 1 );
nroots = mp( length_eta + length_theta + 2 );

br = [ [       q*0;        ones(size(p));   0;   1 ], ...
    [          q*0;         exp(-eta.*y);   0;   0 ], ...
    [    -1./theta;                  p*0;  -1;   0 ], ...
    [          q*0;   eta./(eta - alpha);   0;   1 ], ...
    [          q*0;               1./eta;   0;   0 ] ];

d0 = [0./theta; eta./eta; 0; 1];
d1 = [1./theta; 1./eta; 1; 0];
d2 = [-2*hitlevel./theta - 2./theta.^2; 2./eta.^2; -2*hitlevel; 0];
d3 = [3*hitlevel.^2./theta + 6*hitlevel./theta.^2 + 6./theta.^3; 6./eta.^3; 3*hitlevel.^2; 0];

onesn1 = ones(nroots, 1);
onesn2 = ones(nroots, 1);
onesn1( length_theta + 1 ) = 0;
onesn2( length_theta + 2 ) = 0;

stepsize = 1/mp( 100000 );
load('QvEtau_MC.mat');

vars = [mp(0), 1./mp( [1e15, 1e14, 1e13, 1e12, 1e11, 1e10, 1e9, 1e8, ...
    1e7, 1e6, 1e5, 1e4, 1e3, 1e2] ) ];

betas = mp( vars );
mus = mp( vars );

keepindex = ~ismember(vars, [1e-15, 1e-14, 1e-13, 1e-12, 1e-11]);

mus = unique([mus]); % ; -mus

beta0 = mp( 0 );
mu0 = mp( 0 );

chvar = mp( 0 );
eval( ['chvar = ', chvar_name, 's;'] );
length_chvar = mp( length(chvar) );

J = mp( zeros(length_chvar, 5));
Etau = mp( zeros(length_chvar, 1) );
Qv = mp( zeros(length_chvar, 1) );

for var_th = 1:length_chvar
    
    beta = mp( beta0 );
    mu = mu0;
    eval( [chvar_name, ' = mp( chvar(', num2str(var_th), ') );'] );
    
    if abs(mu)<eps
        mu = mp( 0 );
    end
    if abs(beta)<eps
        beta = mp( 0 );
    end
    
    mubar = mu + lambda * (p'*(1./eta) - q'*(1./theta));
    
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
        A(:, i) = h( Lundbergroots(i), theta, eta, hitlevel );
    end
    
    if beta~=0
        C = A\br;
        J(var_th, :) = exp( - ( hitlevel - x0 ) .* Lundbergroots' ) * C;
        Qv(var_th, :) = ( x0 - hitlevel*J(var_th,1) + J(var_th,3) - J(var_th,5) )./beta ...
            + mubar./beta.^2.*( 1 - J(var_th,1) );
    else
        [ D00, D01, D02, D03, D10, D20, D12, ...
            D10ik, D20ik, D01ik, D02ik, D03ik, D12ik, D20d1k ] = DA( A, d1, d2, d3, br, nroots, length_theta );
        
        ExpRoots_x0 = exp( - ( hitlevel - x0 ) * Lundbergroots' );
        ExpRoots_x0_n1 = ExpRoots_x0 .* onesn1';
        ExpRoots_x0_n2 = ExpRoots_x0 .* onesn2';
        
        if mubar==0
            J(var_th,:) = - ( hitlevel - x0 ) * D10ik( length_theta + 1, : )/D10 ...
                + ExpRoots_x0_n1 * D10ik/D10;
            Etau(var_th,:) = -1/2./zeta .* ( ( hitlevel - x0 ).^2 + ( hitlevel - x0 ) * D20/D10 ...
                + ExpRoots_x0_n1 * D20d1k/D10 );
        elseif mubar~=0
            for i = nroots:-1:1
                A1(:, i) = h1( Lundbergroots(i), mu, sigma, lambda, eta, p, theta, q, hitlevel );
                delta1(i,:) = 1/G1( Lundbergroots(i), mu, sigma, lambda, eta, p, theta, q );
            end
            C = A\br;
            C1 = A\(-A1*C(:,1));
            J(var_th, :) = ExpRoots_x0 * C;
            Etau(var_th,:) = ExpRoots_x0 * ( ( hitlevel - x0 ) .* delta1 .* C(:,1) - C1 );
        end
        
        if mubar>0
            Qv(var_th, :) = 1./mubar/2.*( hitlevel.^2 - x0.^2 + ( 2*hitlevel*D01 + D02 )/D00 ...
                - ExpRoots_x0_n2 * ( 2*hitlevel.*D01ik(:,1) + D02ik(:,1) )/D00 ) ...
                + zeta./mubar.^2.*( x0 - hitlevel + J(var_th,3) - J(var_th,5) );
        elseif mubar==0
            Qv(var_th, :) = 1./zeta/2.*( varpi./zeta.*( ( hitlevel - x0 ).^2 + ( hitlevel - x0 ).*D02/D01 ...
                - ExpRoots_x0_n2 * ( D02ik(:,3) - D02ik(:,5) )/D01 ) ...
                - ( hitlevel - x0 ).*D03./D01/3 - hitlevel.*( hitlevel - x0 + 1 ).*D02/D01 ...
                - ( hitlevel - x0 ).^2.*( x0 + 2*hitlevel )/3 - hitlevel.*D12/D01 ...
                + ExpRoots_x0_n2 * ( 3.*hitlevel.*D12ik(:,1) + D03ik(:,3) - D03ik(:,5) )./D01/3 );
        elseif mubar<0
            Qv(var_th, :) = 1./mubar/2.*( hitlevel.^2 - x0.^2 + ( 2*hitlevel*D10 + D20 )/D00 ...
                - ExpRoots_x0_n1 * ( 2*hitlevel.*D10ik(:,1) + D20ik(:,1) )/D00 ) ...
                + zeta./mubar.^2.*( x0 - hitlevel + J(var_th,3) - J(var_th,5) );
        end
    end
    
    if strcmp( chvar_name, 'beta' )
        fprintf( [chvar_name, ' = %+0.6e, Qv = %0.6e, DiffRoots = %0.6e \r'], ...
            [beta, Qv(var_th,:), Lundbergroots(length_theta + 2) - Lundbergroots(length_theta + 1) ] );
        xlabel_str = '$\beta$';
    else
        fprintf( [chvar_name, ' = %+0.6e, Etau = %0.6e, DiffRoots = %0.6e \r'], ...
            [mu, Etau(var_th,:), Lundbergroots(length_theta + 2) - Lundbergroots(length_theta + 1) ] );
        xlabel_str = '$\bar\mu=\mu$';
    end
end

J2 = exp( - HitTime * diag( sort( [betas, betas] ) ) ) .* ( Overshot > y );

if strcmp( chvar_name, 'beta' )
    J1 = J1(:, 2:2:end);
    J2 = J2(:, 2:2:end);
    J3 = J3(:, 2:2:end);
    J4 = J4(:, 2:2:end);
    J5 = J5(:, 2:2:end);
    QvMC = QvMC2(:, 2:2:end);
    Tau = HitTime(:, 2:2:end);
else
    J1 = J1(:, 1:2:end);
    J2 = J2(:, 1:2:end);
    J3 = J3(:, 1:2:end);
    J4 = J4(:, 1:2:end);
    J5 = J5(:, 1:2:end);
    QvMC = QvMC2(:, 1:2:end);
    Tau = HitTime(:, 1:2:end);
end


CIJ1 = [ nanmean(J1) - 1.96*std(J1)/100; nanmean(J1) + 1.96*std(J1)/100 ]';
CIJ2 = [ nanmean(J2) - 1.96*std(J2)/100; nanmean(J2) + 1.96*std(J2)/100 ]';
CIJ3 = [ nanmean(J3) - 1.96*(std(J3)+sqrt(stepsize))/100; nanmean(J3) + 1.96*(std(J3)+sqrt(stepsize))/100 ]';
CIJ4 = [ nanmean(J4) - 1.96*std(J4)/100; nanmean(J4) + 1.96*std(J4)/100 ]';
CIJ5 = [ nanmean(J5) - 1.96*std(J5)/100; nanmean(J5) + 1.96*std(J5)/100 ]';
CIQv = [ nanmean(QvMC) - 1.96*(std(QvMC)+sqrt(stepsize))/100; nanmean(QvMC) + 1.96*(std(QvMC)+sqrt(stepsize))/100 ]';
CITau = [ nanmean(Tau) - 1.96*(std(Tau)+sqrt(stepsize))/100; nanmean(Tau) + 1.96*(std(Tau)+sqrt(stepsize))/100 ]';

toc

chvar(1) = 1e-15;

Position = [500, 500, 450, 220];

figure(1);
chvartemp = chvar(keepindex);
Jtemp = J(keepindex,1);
CItemp = CIJ1(keepindex,:);
outindex = ( Jtemp < CItemp(:,1) ) + ( Jtemp > CItemp(:,2) ) == 1;
plot( chvartemp, Jtemp, '*-', chvartemp, CItemp, 'm:', 'Linewidth', 2 ); xlabel( xlabel_str, 'Interpreter', 'latex' ); title( '$J_1$', 'Interpreter', 'latex' );
hold on;
plot( chvartemp(outindex), Jtemp(outindex), 'rs', 'Linewidth', 2 );
set( gca, 'xscale', 'log', 'FontSize', 12 );
xticklabel = get( gca, 'xticklabel' );
xticklabel(1) = {'0'};
set( gca, 'xticklabel', xticklabel );
set(gcf, 'Position', Position);
set(gcf,'paperpositionmode','auto');
print(gcf, '-painters', '-depsc','-loose', [OutPutPathFigure,  'J1By_', chvar_name, '_mp.eps'] );

figure(2);
Jtemp = J(keepindex,2);
CItemp = CIJ2(keepindex,:);
outindex = ( Jtemp < CItemp(:,1) ) + ( Jtemp > CItemp(:,2) ) == 1;
plot( chvartemp, Jtemp, '*-', chvartemp, CItemp, 'm:', 'Linewidth', 2 ); xlabel( xlabel_str, 'Interpreter', 'latex' ); title( '$J_2$', 'Interpreter', 'latex' );
hold on;
plot( chvartemp(outindex), Jtemp(outindex), 'rs', 'Linewidth', 2 );
set( gca, 'xscale', 'log', 'FontSize', 12 );
xticklabel = get( gca, 'xticklabel' );
xticklabel(1) = {'0'};
set( gca, 'xticklabel', xticklabel );
set(gcf, 'Position', Position);
set(gcf,'paperpositionmode','auto');
print(gcf, '-painters', '-depsc','-loose', [OutPutPathFigure,  'J2By_', chvar_name, '_mp.eps'] );
close all;

figure(2);
Jtemp = J(keepindex,2);
CItemp = CIJ2(keepindex,:);
outindex = ( Jtemp < CItemp(:,1) ) + ( Jtemp > CItemp(:,2) ) == 1;
plot( chvartemp, 1-Jtemp, '*-', chvartemp, 1-CItemp, 'm:', 'Linewidth', 2 ); xlabel( xlabel_str, 'Interpreter', 'latex' ); title( '\bf{Q (ii)}', 'Interpreter', 'latex' );
hold on;
plot( chvartemp(outindex), 1-Jtemp(outindex), 'rs', 'Linewidth', 2 );
set( gca, 'xscale', 'log', 'FontSize', 12 );
xticklabel = get( gca, 'xticklabel' );
xticklabel(1) = {'0'};
set( gca, 'xticklabel', xticklabel );
set(gcf, 'Position', Position);
set(gcf,'paperpositionmode','auto');
print(gcf, '-painters', '-depsc','-loose', [OutPutPathFigure,  'Q2By_', chvar_name, '_mp.eps'] );

figure(3);
Jtemp = J(keepindex,3);
CItemp = CIJ3(keepindex,:);
outindex = ( Jtemp < CItemp(:,1) ) + ( Jtemp > CItemp(:,2) ) == 1;
plot( chvartemp, Jtemp, '*-', chvartemp, CItemp, 'm:', 'Linewidth', 2 ); xlabel( xlabel_str, 'Interpreter', 'latex' ); title( '\bf{Q (iv)}', 'Interpreter', 'latex' );
hold on;
plot( chvartemp(outindex), Jtemp(outindex), 'rs', 'Linewidth', 2 );
set( gca, 'xscale', 'log', 'FontSize', 12 );
xticklabel = get( gca, 'xticklabel' );
xticklabel(1) = {'0'};
set( gca, 'xticklabel', xticklabel );
set(gcf, 'Position', Position);
set(gcf,'paperpositionmode','auto');
print(gcf, '-painters', '-depsc','-loose', [OutPutPathFigure,  'Q4By_', chvar_name, '_mp.eps'] );

figure(4);
Jtemp = J(keepindex,4);
CItemp = CIJ4(keepindex,:);
outindex = ( Jtemp < CItemp(:,1) ) + ( Jtemp > CItemp(:,2) ) == 1;
plot( chvartemp, Jtemp, '*-', chvartemp, CItemp, 'm:', 'Linewidth', 2 ); xlabel( xlabel_str, 'Interpreter', 'latex' ); title( '$J_4$', 'Interpreter', 'latex' );
hold on;
plot( chvartemp(outindex), Jtemp(outindex), 'rs', 'Linewidth', 2 );
set( gca, 'xscale', 'log', 'FontSize', 12 );
xticklabel = get( gca, 'xticklabel' );
xticklabel(1) = {'0'};
set( gca, 'xticklabel', xticklabel );
set(gcf, 'Position', Position);
set(gcf,'paperpositionmode','auto');
print(gcf, '-painters', '-depsc','-loose', [OutPutPathFigure,  'J4By_', chvar_name, '_mp.eps'] );

figure(5);
Jtemp = J(keepindex,5);
CItemp = CIJ5(keepindex,:);
outindex = ( Jtemp < CItemp(:,1) ) + ( Jtemp > CItemp(:,2) ) == 1;
plot( chvartemp, Jtemp, '*-', chvartemp, CItemp, 'm:', 'Linewidth', 2 ); xlabel( xlabel_str, 'Interpreter', 'latex' ); title( '$J_5$', 'Interpreter', 'latex' );
hold on;
plot( chvartemp(outindex), Jtemp(outindex), 'rs', 'Linewidth', 2 );
set( gca, 'xscale', 'log', 'FontSize', 12 );
xticklabel = get( gca, 'xticklabel' );
xticklabel(1) = {'0'};
set( gca, 'xticklabel', xticklabel );
set(gcf, 'Position', Position);
set(gcf,'paperpositionmode','auto');
print(gcf, '-painters', '-depsc','-loose', [OutPutPathFigure,  'J5By_', chvar_name, '_mp.eps'] );

figure(6);
Jtemp = Qv(keepindex);
CItemp = CIQv(keepindex,:);
outindex = ( Jtemp < CItemp(:,1) ) + ( Jtemp > CItemp(:,2) ) == 1;
plot( chvartemp, Jtemp, '*-', chvartemp, CItemp, 'm:', 'Linewidth', 2 ); xlabel( xlabel_str, 'Interpreter', 'latex' ); title( '\bf{Q (v)}', 'Interpreter', 'latex' );
hold on;
plot( chvartemp(outindex), Jtemp(outindex), 'rs', 'Linewidth', 2 );
set( gca, 'xscale', 'log', 'FontSize', 12 );
xticklabel = get( gca, 'xticklabel' );
xticklabel(1) = {'0'};
set( gca, 'xticklabel', xticklabel );
set(gcf, 'Position', Position);
set(gcf,'paperpositionmode','auto');
print(gcf, '-painters', '-depsc','-loose', [OutPutPathFigure,  'Q5By_', chvar_name, '_mp.eps'] );

if strcmp( chvar_name, 'mu' ) && beta0 == 0
    figure(7);
    Jtemp = Etau(keepindex);
    CItemp = CITau(keepindex,:) + ( CITau(keepindex,:) - mean( CITau(keepindex,:), 2 ) );
    outindex = ( Jtemp < CItemp(:,1) ) + ( Jtemp > CItemp(:,2) ) == 1;
    plot( chvartemp, Jtemp, '*-', chvartemp, CItemp, 'm:', 'Linewidth', 2 ); xlabel( xlabel_str, 'Interpreter', 'latex' );
    title( '\bf{Q (iii)}', 'Interpreter', 'latex' ); % $E_x[\tau_b]$
    hold on;
    plot( chvartemp(outindex), Jtemp(outindex), 'rs', 'Linewidth', 2 );
    set( gca, 'xscale', 'log', 'FontSize', 12 );
    xticklabel = get( gca, 'xticklabel' );
    xticklabel(1) = {'0'};
    set( gca, 'xticklabel', xticklabel );
    set(gcf, 'Position', Position);
    set(gcf,'paperpositionmode','auto');
    print(gcf, '-painters', '-depsc','-loose', [OutPutPathFigure,  'Q3By_', chvar_name, '_mp.eps'] );
end

return;

end


%%%%%%%%%%%%%%%%%%%% Function h ( Columns of Matrix A(beta) ) %%%%%%%%%%%%%%%%%%%%%%%%
function [ hx ] = h( x, theta, eta, hitlevel )
hx = [ [ x ./ ( theta + x ) .* exp( - x .* hitlevel ) ]; ...
    [                              eta ./ ( eta - x ) ]; ...
    [                     x .* exp( - x .* hitlevel ) ]; ...
    [                                               1 ];     ]; %#ok<*NBRAK>
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
function [] = RHEP_MC(  )

clc;
close all;

x0 = 0;
mu = 0;
sigma = 0.3;
lambda = 5;
hitlevel = 0.1; % hitting level
y = 0.05; % overshot level;

eta = linspace( 20, 40, 10 )'; % exp parameters for upward jumps
p = ones(size(eta))/length( eta )/2; % probability for upward jumps
theta = eta; % exp parameters for downward jumps
q = p; % probability for downward jumps
alpha = min(eta)/2; % For J_4 ( Laplace Transform Parameter )

quantile_hep = Find_quantile_hep( eta, p, theta, q );

t_all = 2;
stepsize = 1/100000;
ts = ( 0:stepsize:t_all )';
Nstep = round( t_all/stepsize );
Nsample = 100;
Nexperiment = 100;

sigmadt = sigma * sqrt(stepsize);

vars = [ 0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2 ];
%return;

ujump = rand( Nstep, Nsample );
ujump = min( 1e6, floor( ujump*1e6 ) + 1 );
ujump = quantile_hep( ujump );

QvMC2 = zeros( Nexperiment*Nsample, length(vars)*2 )/0;
HitTime = zeros( Nexperiment*Nsample, length(vars)*2 )/0;
Overshot = zeros( Nexperiment*Nsample, length(vars)*2 )/0;
J1 = zeros( Nexperiment*Nsample, length(vars)*2 )/0;
J2 = zeros( Nexperiment*Nsample, length(vars)*2 )/0;
J3 = zeros( Nexperiment*Nsample, length(vars)*2 )/0;
J4 = zeros( Nexperiment*Nsample, length(vars)*2 )/0;
J5 = zeros( Nexperiment*Nsample, length(vars)*2 )/0;

TrueValue = 0.003491562;
t0 = now;
for i = 1:Nexperiment
    RandNormal = randn( Nstep, Nsample );
    RandUnit = rand( Nstep, Nsample );
    ulambda = RandUnit <= lambda * stepsize;
    
    for var_th = 1:length(vars)
        for mb_th = 1:2
            if mb_th == 1
                mu = vars(var_th);
                beta = 0;
            else
                mu = 0;
                beta = vars(var_th);
            end
            
            ExpBetaTs = exp( - beta * ts );
            wn = sigmadt * RandNormal + mu * stepsize;
            
            Xt = cumsum( wn + ulambda.*ujump );
            Xt = [Xt(1,:)*0; Xt]; %#ok<AGROW>
            Lt = - min( 0, cummin(Xt) );
            XTt = Xt + Lt;
            
            %%% In light of Broadie and Glasserman (MF, 1997), Broadie,
            %%% Glasserman, Kou (SF, 1999), we need to correct the barrier
            %%% by a factor of exp( 0.5826 * sigmadt )
            Adj_factor = exp( - 0.5826 * sigmadt );
            hitlevel_Adj = hitlevel * Adj_factor;
            
            % return;
            PreTauIndex = cumsum( XTt > hitlevel_Adj ) == 0;
            XTtPreTau = XTt.*PreTauIndex;
            CumMaxXTt = cummax( XTt );
            
            HitTime_temp = sum( CumMaxXTt < hitlevel_Adj ) * stepsize;
            HitTime( (i-1)*Nsample+1:i*Nsample, (var_th-1)*2 + mb_th ) = HitTime_temp;
            
            CumMaxXTt_temp = CumMaxXTt;
            CumMaxXTt_temp(PreTauIndex) = nan;
            Overshot_temp = nanmin( CumMaxXTt_temp ) - hitlevel_Adj;
            
            [PostTauValue, row_indices] = nanmin( CumMaxXTt_temp );
            % linearInd_CumMaxXTt_temp = sub2ind( size(CumMaxXTt_temp), max(1, row_indices - 1), 1:size(CumMaxXTt_temp, 2));
            % PreTauValue = CumMaxXTt(linearInd_CumMaxXTt_temp);
            linearInd_ulambda = sub2ind( size(ulambda), max(1, row_indices - 1), 1:size(CumMaxXTt_temp, 2));
            JumpOverIndex = ulambda(linearInd_ulambda);
            % fprintf( '%0.6f,  %0.6f,  %3d,  %8d,  %3d\r', [PreTauValue', PostTauValue', JumpOverIndex', row_indices' - 1, (1:size(CumMaxXTt_temp, 2))']' );
            % return;
            
            Overshot_temp = Overshot_temp.*JumpOverIndex;
            Overshot( (i-1)*Nsample+1:i*Nsample, (var_th-1)*2 + mb_th ) = Overshot_temp;
            
            % QvMC1(i,:) = ExpBetaTs' * XTtPreTau * stepsize;
            %%%% next is similar to trapz rule
            QvMC2( (i-1)*Nsample+1:i*Nsample, (var_th-1)*2 + mb_th ) ...
                = ( ExpBetaTs(1:end-1)' * XTtPreTau(1:end-1, :) + ExpBetaTs(2:end)' * XTtPreTau(2:end, :) )/2 * stepsize ...
                - hitlevel_Adj * stepsize/2;
            
            J1( (i-1)*Nsample+1:i*Nsample, (var_th-1)*2 + mb_th ) = exp( - HitTime_temp * beta );
            J2( (i-1)*Nsample+1:i*Nsample, (var_th-1)*2 + mb_th ) = exp( - HitTime_temp * beta ) .* ( Overshot_temp > y );
            J3( (i-1)*Nsample+1:i*Nsample, (var_th-1)*2 + mb_th ) ...
                = ExpBetaTs(1:end-1)' * ( diff( Lt ) .* PreTauIndex(1:end-1,:) );
            J4( (i-1)*Nsample+1:i*Nsample, (var_th-1)*2 + mb_th ) = exp( - HitTime_temp * beta + alpha * Overshot_temp );
            J5( (i-1)*Nsample+1:i*Nsample, (var_th-1)*2 + mb_th ) = exp( - HitTime_temp * beta ) .* Overshot_temp;
            
            %             [nanmean(J1(:)), nanstd(J1(:))/sqrt(Nsample)]
            %             [nanmean(J2(:)), nanstd(J2(:))/sqrt(Nsample)]
            %             [nanmean(J3(:)), nanstd(J3(:))/sqrt(Nsample)]
            %             [nanmean(J4(:)), nanstd(J4(:))/sqrt(Nsample)]
            %             [nanmean(J5(:)), nanstd(J5(:))/sqrt(Nsample)]
            %             [mean(Overshot), std(Overshot)/sqrt(Nsample)]
            %             return;
            
            data1 = QvMC2( 1:i*Nsample, (var_th-1)*2 + mb_th );
            data1 = data1(:);
            data2 = HitTime( 1:i*Nsample, (var_th-1)*2 + mb_th );
            data2 = data2(:);
            
            MeanQvMC = nanmean(data1);
            StdQvMC = nanstd(data1)/sqrt(i*Nsample);
            fprintf( 'N = (%03d, %02d, %02d), TrueValue = %0.6f, MC = (%0.6f, %0.6f), MeanQvMC = %0.6f, StdQvMC = %0.6f, MaxTau = %0.6f, TimeUsed = %03.2f \r', ...
                [i, var_th, mb_th, TrueValue, MeanQvMC - 1.96 * StdQvMC, nanmean(data1) + 1.96 * StdQvMC, ...
                MeanQvMC, StdQvMC, nanmax( data2 ), ( now - t0 ) * 86400] );
            
        end
    end
    
end

save( 'QvEtau_MC.mat', 'HitTime', 'QvMC2', ...
    'Overshot', 'J1', 'J2', 'J3', 'J4', 'J5', 'vars' );

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
    if mod(i, 1000) == 1
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
