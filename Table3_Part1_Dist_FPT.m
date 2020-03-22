function [] = Table3_Part1()
clc;

mu = 0.1;
sig = 0.3;
eta = [20; 30];
theta = [30; 40];
p = [0.2; 0.3];
q = [0.2; 0.3];
x = 0;
A0 = 20;
n0 = 10;

t = 0.5;
ld = 2;
b = 0.3;

tic;
LIV = Euler_cdf_fps_ref_hem_general(mu,sig,ld,eta,theta,p,q,x,b,t,A0,n0);
cputime = toc;

fprintf('t = %0.2f, lambda = %0.2f, b = %0.2f, LIV = %0.5f, CUP Time = %0.4f sec.\r', ...
    [t, ld, b, LIV, cputime]);

return;

end

%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Below are subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
function f=Euler_cdf_fps_ref_hem_general(mu,sig,ld,eta,theta,p,q,x,b,t,A1,n1)
%This program is to compute the distribution of the first passage time of
%the reflected hyper-exponential jump diffusion process (RHEP) via the Euler
%inversion algorithm, i.e., P(\tau_b\leq t)

%Input:
%     mu: drift
%     sig: volatility
%     ld: jump rate
%     eta: exponential rates of the upward hyper-exponential jump
%     theta: exponential rates of the downward hyper-exponential jump 
%     p: probability weights of the upward hyper-exponential jump
%     q: probability weights of the downward hyper-exponential jump
%     x: starting point of the RHEP
%     b: the barrier
%     t: the time horizon
%     A1: a parameter of the Euler inversion algorithm  
%     n1: a parameter of the Euler inversion algorithm
%
%Output:
%     f: the value of P(\tau_b\leq t)
   
m1=n1+15;
coef=exp(A1/2)/(2*t);
g=zeros(1,m1+n1+1);

for j=1:(m1+n1+1)
    temp=Lap_cdf_fpt_g_reflected_hem_any(mu,sig,ld,eta,theta,p,q,x,b,(A1-2*(j-1)*pi*i)/(2*t));
    g(j)=real(temp);
end

f = Euler(n1, m1, A1, g);
f = 2*coef*f-coef*g(1);
end

%%%%%%%%%%%%%%%
function f=Lap_cdf_fpt_g_reflected_hem_any(mu,sig,ld,eta,theta,p,q,x,b,s)

%This program is to calculate the Laplace transform of the cdf of the first 
%passage time of the reflected hyper-exponential jump diffusion process.

%Input:
%     mu: drift
%     sig: volatility
%     ld: jump rate
%     eta: exponential rates of the upward hyper-exponential jump
%     theta: exponential rates of the downward hyper-exponential jump 
%     p: probability weights of the upward hyper-exponential jump
%     q: probability weights of the downward hyper-exponential jump
%     x: starting point of the RHEP
%     b: the barrier
%     s: the point at which the Laplace transform is computed
%
%Output:
%     f: the Laplace transform (at the point s) of the cdf

f=Lap_pdf_fpt_g_reflected_hem_any(mu,sig,ld,eta,theta,p,q,x,b,s)/s;
end

%%%%%%%%%%%%%%%
function f = Euler(n, m, A, TransformedFValue)
%Euler_Sum returns an approximate sum of a series like 
%      \sum_i=0^infty (-1)^i * TransformedFValue(i+1)
%using Euler_Sum method.
%Note that TransformedFvalue(i) is the (i-1)th value of the original series.
% TransformedFValue is a vector with m+n+1 elements:
% TransformedFValue=(g(A/2*t), g((A+2*pi*i)/(2*t)), g((A+2*2*pi*i)/(2*t)), ..., g((A+2*(m+n)*pi*i)/(2*t))',
% where g is the Laplace transformation of f. 

S=zeros(m+1,1);
S(1)=TransformedFValue(1);

sign=-1;
for k=1:n
    S(1)=S(1)+sign*TransformedFValue(k+1);
    sign=sign*(-1);
end
sign = (-1)^(n+1);
for j=1:m
    S(j+1)=S(j)+sign*TransformedFValue(n+j+1);
    sign=sign*(-1);
end
coefs=binopdf((0:m)',m,0.5);
f=sum(S.*coefs);
end

%%%%%%%%%%%%%%%
function f=Lap_pdf_fpt_g_reflected_hem_any(mu,sig,ld,eta,theta,p,q,x,b,s)

%This program is to calculate the Laplace transform of the pdf of the first 
%passage time of the reflected hyper-exponential jump diffusion process.

%Input:
%     mu: drift
%     sig: volatility
%     ld: jump rate
%     eta: exponential rates of the upward hyper-exponential jump
%     theta: exponential rates of the downward hyper-exponential jump 
%     p: probability weights of the upward hyper-exponential jump
%     q: probability weights of the downward hyper-exponential jump
%     x: starting point of the RHEP
%     b: the barrier
%     s: the point at which the Laplace transform is computed
%
%Output:
%     f: the Laplace transform (at the point s) of the pdf

pu=sum(p);
qd=sum(q);

delta=roots_Gx_g_reflected_hem(mu,sig,ld,eta,theta,pu,p/pu,qd,q/qd,s);
m=length(eta);n=length(theta);


A=ones(n+m+2,n+m+2);
for j=1:(n+m+2)
    A(:,j)=h_reflected_hem(eta, theta, b, delta(j));
end

b1=ones(n+m+2,1);
b1(1:n)=0;
b1(n+m+1)=0;
C1=A\b1;
c=C1.*exp(-delta*(b-x));
f=sum(c);
end

%%%%%%%%%%%%%%%
function delta=roots_Gx_g_reflected_hem(mu,sig,ld,eta,theta,pu,p,qd,q,s)
%This program is to calculate all the roots of function G(x)=s
%with any dimension of eta and theta

%Input:
%     mu: drift
%     sig: volatility
%     ld: jump rate
%     eta: exponential rates of the upward hyper-exponential jump
%     theta: exponential rates of the downward hyper-exponential jump 
%     pu: upward jump probability 
%     p: probability weights of the upward hyper-exponential jump / pu
%     qd: downward jump probability
%     q: probability weights of the downward hyper-exponential jump / qd
%     s: the "s" involved in the function G(x)=s
%
%Output:
%     delta: all the roots of function G(x)=s

newsig=sig*sig/2;

m=length(eta);n=length(theta);
eta=reshape(eta,1,m);theta=reshape(theta,1,n);
p=reshape(p,1,m);q=reshape(q,1,n);
k=m+n;
pi=[-eta theta];
r=[pu*p qd*q];

a=zeros(k,k+3);

a(1,1)=ld*r(1)*pi(1)-pi(1)*(ld+s);
a(1,2)=mu*pi(1)-ld-s;
a(1,3)=newsig*pi(1)+mu;
a(1,4)=newsig;

if k>1
    b=zeros(k-1,k);
    b(1,1)=pi(1); 
    b(1,2)=1;
    if k>2
        for l=2:(k-1)
            b(l,1)=b(l-1,1)*pi(l);
            b(l,l+1)=b(l-1,l);
            for h=2:l
                b(l,h)=pi(l)*b(l-1,h)+b(l-1,h-1);
            end
        end
    end
    for j=2:k
        a(j,1)=a(j-1,1)*pi(j)+b(j-1,1)*ld*r(j)*pi(j);
        a(j,j+3)=a(j-1,j+2);
        a(j,j+2)=a(j-1,j+2)*pi(j)+a(j-1,j+1);
        a(j,j+1)=a(j-1,j+1)*pi(j)+a(j-1,j);
        for i=2:j
            a(j,i)=a(j-1,i-1)+a(j-1,i)*pi(j)+b(j-1,i)*ld*r(j)*pi(j);
        end
    end
end

c=a(k,((k+3):(-1):1));
allroot=roots(c);
[temp index]=sort(real(allroot));
delta=allroot(index);
end

%%%%%%%%%%%%%%%
function f=h_reflected_hem(eta,theta,b,x)
%This program is to compute the functon h(.) for the purpose of forming 
%the matrix A(beta) under the reflected hyper-exponential jump diffusion 
%process

%Input:
%     eta: exponential rates of the upward hyper-exponential jump
%     theta: exponential rates of the downward hyper-exponential jump 
%     b: the barrier
%     x: the point at which the function h(.) is evaluated.
%
%Output:
%     f: value of the function h(x)

m=length(eta);
n=length(theta);

f=zeros(m+n+2,1);

f(1:n)=x*exp(-x*b)./(theta+x);
f((n+1):(n+m))=eta./(eta-x);
f(n+m+1)=x*exp(-x*b);
f(n+m+2)=1;
end

%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%% The End %%%%%%%%%%%%%%%