function l=LoglikelihoodRiemannDerivatives(T,spiketrain,spikebins,lags,f,g,tau,tau2,M,type)

%Type=1 corresponds to DR1
%Type=2 corresponds to DR2
%Define bins
points=linspace(0,T,M+1);
delta=T/M;
%Define evaluation grid
if(tau2>0)
h=@(x)1/tau2*(x-tau);
r=@(x)max(0,min(h(x),1));
else
%if tau2=0 then r(t)=1 for t>tau (used for the renewal
%process case)
r=@(x)x>tau;
end

%evaluate sum_i Delta_i g(t_i) - f(t_i)r(t_i-t_N(t_i))*delta
%with Delta_i the number of spikes in bin i
%For the log likelihood  evaluations, g(x)=log(f(x))
%For the rest of evaluations f and g are the derivatives of
%the intensity and log intensity, respectively.

vec1=arrayfun(g,points(1:end-1));
vec1(isnan(vec1))=0;
vec1(isinf(vec1))=0;
vec2=arrayfun(f,points(1:end-1)).*arrayfun(r,lags)*delta;
l=dot(spikebins,vec1)-sum(vec2);

if(type==1)
 
    return
else
   l=l+dot(spikebins,vec2)/2; return
end
    