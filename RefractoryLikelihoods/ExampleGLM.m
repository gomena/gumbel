disp('This Example simulates a spike train and then does ')
disp('parameter estimation using Discrete Riemann 1 (DR1)') 
disp('and Gauss-Lobatto (GL) approximations' )

%In this case the intensity function is given by
%f(x)=exp(K1*sin(2*pi*lx)+K2), K=[K1,K2] is the parameter to be inferred
%the conditional intensity function is
%f2(x)=f(x)*r(x-t_N(t)-)
%r(t)=0 for t in [0,tau]
%r(t) increases from 0 to 1 between [tau,tau+tau2]
%r(t)=1 for t>tau+tau2

%First load the pre-computed quadrature nodes and weigths

%Define functions and parameters
K1=3; K2=2;
K=[K1 K2]';
tau=0.002;
tau2=0.01;
l=2;
T=40;

disp(['The intensity function is f(t)=exp(K1*sin(2pi*l*t)+K2)'])
disp(['With l=' num2str(l) ',K1=' num2str(K1) '(Gain) and K2=' num2str(K2) '(baseline)'])
disp(['K=[K1 K2]=[' [num2str(K')] '] is the parameter vector to be estimated.']) 
disp(['Now Simulating a Spike Train for T = ' num2str(T) ' seconds'])


delta=0.0001;
g=@(x)1/tau2*(x-tau);
r=@(x)max(0,min(g(x),1));
f=@(x)exp(sin(l*x*2*pi)*K1+K2);
pause(0.5)
disp(['To account for refractoriness, the actual intensity'])
disp(['is f2(t)=f(t)*r(t-N(t)-) , r(t)=0 in [0,tau], r(t) goes from'])
disp(['0 to 1 in [tau,tau+tau2] and r(t)=1 after tau+tau2'])
disp(['Here, tau=' num2str(tau) ' and tau2=' num2str(tau2)])
%Simulate spiketrain using the time-change theorem 
%(see for example Citi et al 2014)
t=cputime;
  spiketrain=[];
  flag=0;
  times=[];
  feval=[] ;
while(flag==0)  
    sum=0;
    sample=-log(unifrnd(0,1));
    i=0;sum=0;
    
    while(sum<sample)
        if(isempty(spiketrain))
    sum=sum+f(delta*i)*delta; 
    i=i+1; 
        else   
            
    sum=sum+f(delta*i+spiketrain(end))*delta*r(delta*i);
    i=i+1;
        end
    end
    if(isempty(spiketrain))
        spiketrain=[spiketrain i*delta];
    else
    spiketrain=[spiketrain spiketrain(end)+i*delta];
    end
    if(max(spiketrain)>T)
        flag=1;
    end
end
spiketrain=spiketrain(find(spiketrain<T-tau));
disp(['It Took ' num2str(cputime-t) ' seconds to simulate the spike train']) 

%The following code will generate two plots
%The first is the intensity between T=0 and T=5;
%The second plot is the conditional
%intensity function during the first 40 spikes (as in IntensityGLM.png)

subplot(1,2,1)
plot([0:0.01:5],arrayfun(@(x)f(x),[0:0.01:5]))
xlim([0 5])
subplot(1,2,2)
spikeplot=[spiketrain(1:6)];
for i=1:length(spikeplot)-1
    hold on 
plot(spikeplot(i+1),f(spikeplot(i+1))*r(spikeplot(i+1)-spikeplot(i)),'ko','MarkerFaceColor','k','MarkerSize',3,'LineWidth',1.5)
points=[spikeplot(i):0.0001:spikeplot(i+1)];
plot(points,f(points).*r(points-points(1)),'black','Linewidth',1);
end

xlim([0 spikeplot(end)])
disp('Two plots were generated. The left plot is')
disp('f(t) between T=0 and T=5')
disp('The right plot is the CIF for the first 6 spikes')

disp('Now we will estimate the parameters maximizing the likelihood')
disp('Using Newton method with backtracking linesearches')



%Evaluations of the likelihood per integral
M=10*length(spiketrain);
disp(['A number of' num2str(M) ' likelihood evaluations will be used'])
disp(['This is ten times the number of spikes']) 
disp(['You can try different values by tweaking the code'])
%Minimum number of evaluations per integral
mineval=1;
p=90;

%Now define again the intensity, this time 
%also as a function of the unknown parameter k
f=@(x,k)exp(sin(l*x*2*pi)*k(1)+k(2));
%Same with the log intensity
logf=@(x,k)k(1)*sin(l*x*2*pi)+k(2);

%Same with derivatives of intensity and log intensity
df1=@(x,k)sin(l*x*2*pi)*exp(k(1)*sin(l*x*2*pi)+k(2));
df2=@(x,k)exp(sin(l*x*2*pi)*k(1)+k(2));

d2f11=@(x,k)sin(l*x*2*pi)^2*exp(k(1)*sin(l*x*2*pi)+k(2));
d2f12=@(x,k)sin(l*x*2*pi)*exp(k(1)*sin(l*x*2*pi)+k(2));
d2f22=@(x,k)exp(k(1)*sin(l*x*2*pi)+k(2));

dlogf1=@(x,k)sin(l*x*2*pi);
dlogf2=@(x,k)1;


d2logf11=@(x,k)0;
d2logf12=@(x,k)0;
d2logf22=@(x,k)0;


shg
disp('Press any key to begin')
pause
tol=1e-4; %Gradient tolerance  for all Newton optimizations
maxiter=6; %maximum number of iteartions for all newton optimizations
k0=[0 0]'; %Initial point for all Newton optimizations

disp('Now doing optimization with Riemann (DR1) approximation')
t=cputime;

%Define grid and bins to evaluate the formula
spikebins=BinSpikes(spiketrain,T,M);
points=linspace(0,T,M+1);
%Find vectors of lags from the last spike to compute
%refractory effects
lags=tLastSpike(points(1:end-1),spiketrain);
%Define functions
type=1;
logfunriemann=@(k)LoglikelihoodRiemannDerivatives(T,spiketrain,spikebins,lags,@(x)f(x,k),@(x)logf(x,k),tau,tau2,M,type);
dlogfunriemann1=@(k)LoglikelihoodRiemannDerivatives(T,spiketrain,spikebins,lags,@(x)df1(x,k),@(x)dlogf1(x,k),tau,tau2,M,type);
dlogfunriemann2=@(k)LoglikelihoodRiemannDerivatives(T,spiketrain,spikebins,lags,@(x)df2(x,k),@(x)dlogf2(x,k),tau,tau2,M,type);
gradlogfunriemann=@(k)[dlogfunriemann1(k) dlogfunriemann2(k)]'; 

d2logfunriemann11=@(k)LoglikelihoodRiemannDerivatives(T,spiketrain,spikebins,lags,@(x)d2f11(x,k),@(x)d2logf11(x,k),tau,tau2,M,type);
d2logfunriemann12=@(k)LoglikelihoodRiemannDerivatives(T,spiketrain,spikebins,lags,@(x)d2f12(x,k),@(x)d2logf12(x,k),tau,tau2,M,type);
d2logfunriemann22=@(k)LoglikelihoodRiemannDerivatives(T,spiketrain,spikebins,lags,@(x)d2f22(x,k),@(x)d2logf22(x,k),tau,tau2,M,type);
Hlogfunriemann=@(k)[[d2logfunriemann11(k) d2logfunriemann12(k)];[d2logfunriemann12(k) d2logfunriemann22(k)]];

%Do optimization using Newton method with Backtrack linesearches
KDR1=NewtonBackTrack(logfunriemann,gradlogfunriemann,Hlogfunriemann,tol,maxiter,k0);


disp(['Optimization finished. Took ' num2str(cputime-t) ' seconds']) 
disp('The vector KDR1 contains the estimated parameters using the DR1 approximation')


disp('Now doing optimization with the Gauss-Lobatto approximation')



t=cputime;
%Find nodes and weights
[node weight neval]=NodesWeightsGL(spiketrain,T,tau,M);
%Find vectors of lags from the last spike to compute
%refractory effects
lagsGL=tLastSpike(node,spiketrain);
%define functions
logfunGL=@(k)LoglikelihoodGLDerivatives(T,spiketrain,lagsGL,@(x)f(x,k),@(x)logf(x,k),tau,tau2,weight,node);
dlogfunGL1=@(k)LoglikelihoodGLDerivatives(T,spiketrain,lagsGL,@(x)df1(x,k),@(x)dlogf1(x,k),tau,tau2,weight,node);
dlogfunGL2=@(k)LoglikelihoodGLDerivatives(T,spiketrain,lagsGL,@(x)df2(x,k),@(x)dlogf2(x,k),tau,tau2,weight,node);
gradlogfunGL=@(k)[dlogfunGL1(k) dlogfunGL2(k)]'; 



d2logfunGL11=@(k)LoglikelihoodGLDerivatives(T,spiketrain,lagsGL,@(x)d2f11(x,k),@(x)d2logf11(x,k),tau,tau2,weight,node);
d2logfunGL12=@(k)LoglikelihoodGLDerivatives(T,spiketrain,lagsGL,@(x)d2f12(x,k),@(x)d2logf12(x,k),tau,tau2,weight,node);
d2logfunGL22=@(k)LoglikelihoodGLDerivatives(T,spiketrain,lagsGL,@(x)d2f22(x,k),@(x)d2logf22(x,k),tau,tau2,weight,node);
HlogfunGL=@(k)[[d2logfunGL11(k) d2logfunGL12(k)];[d2logfunGL12(k) d2logfunGL22(k)]];




%Do optimization using Newton method with Backtrack linesearches
KGL=NewtonBackTrack(logfunGL,gradlogfunGL,HlogfunGL,tol,maxiter,k0);

disp(['Optimization finished. Took ' num2str(cputime-t) ' seconds']) 

disp(['KDR1 = [' [num2str(KDR1')] '] is the estimate obtained using DR1']);
disp(['KGL  = [' [num2str(KGL')]  '] is the estimate obtained using Gauss-Lobatto']);

disp('Recall that as M goes to infinity both approximations converge')
disp(['To the ML estimator of K = [' [num2str(K')] '] '])
disp('which is different from K (but converges to K as T goes to infinity)')
close(gcf)