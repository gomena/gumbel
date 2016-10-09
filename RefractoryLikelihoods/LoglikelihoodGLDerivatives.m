function [y ints]=LoglikelihoodGLDerivatives(T,spiketrain,lags,f,g,tau,tau2,w,x)


%First define refractory Kernell
%r(t)=0 for t in [0,tau]
%r(t) increases from 0 to 1 between [tau,tau+tau2]
%r(t)=1 for t>tau+tau2

if(tau2>0)
h=@(x)1/tau2*(x-tau);
r=@(x)max(0,min(h(x),1));
else
%if tau2=0 then r(t)=1 for t>tau (used for the renewal
%process case)
r=@(x)x>tau;
end
%Compute sum_i g(t_i) - sum_k f(t_k)w_k
%where t_i are spiking times
%t_k the quadrature nodes
%and w_k the quadrature weights
%(recall that t_k can be chosen so that t_i are a subset of
%t_k)
%For the log likelihood  evaluations, g(x)=log(f(x))
%For the rest of evaluations f and g are the derivatives of
%the intensity and log intensity, respectively.

y=sum(arrayfun(g,spiketrain))-dot(w,arrayfun(f,x).*arrayfun(r,lags));