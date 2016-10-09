%In this example we compute the log likelihood of an Inverse Gaussian
%distribution

%Load nodes and weights

disp('In this example a Renewal Process with an Inverse Gaussian')
disp('ISI-distribution will be simulated to compare (against the true value)')
disp('the four approximations considered in the manuscript:')
disp('two based on binary time series, discrete Riemann 1 and 2 (DR1, DR2)')
disp('and two based on the continuous time process; continuous Trapezoidal (CT)')
disp('and Gauss-Lobatto (GL)')


%Define parameters lambda, mu of the inv gauss distribution 

mu=1/10;
lambda=1;
%Define refractory period
tau=0.002;
%Define recording time, 100s
T=100;
disp(['The simulation time is T=' num2str(T) '(s) and the parameters'])
disp(['of the ISI dist. are mu=' num2str(mu) ',lambda=' num2str(lambda) ])
disp(['and absolute refractory period tau=' num2str(tau)])
disp('Press any key to simulate spike train')
pause

%Define the non-translated ISI distribution and draw samples
PD = ProbDistUnivParam('inversegaussian',[mu lambda]);
%draw samples
y=random(PD,10000,1);


%Define the corresponding CDF,PDF, and Intensity function
cdfinvgauss=@(x)cdf(PD,x);
pdfinvgauss=@(x)pdf(PD,x);
intinvgauss=@(x)pdf(PD,x)./(1-cdf(PD,x));
%Define the true ISI distribution and compute the CDF
%PDF and intensity functions (the actual
%ISI distribution is a tau-translated version of the
%inverse gaussian
[pdftrans,cdftrans,inttrans]=TransFunctions(pdfinvgauss,cdfinvgauss,intinvgauss,tau);


%Add the refractory period to the samples 
y=y+tau;

%Compute the spike train. 
%Assume, for convenience, the first spike occurs right before t=0

spiketrain(1)=y(1);
for i=2:length(y)
spiketrain(i)=spiketrain(i-1)+y(i);
end

spiketrain=spiketrain(find(spiketrain<T-tau));

int=@(x)RenewalIntensity(x,inttrans,spiketrain); %intensity function
tau2=0; %(no refractory function r, see LoglikelihoodGLDerivatives.m for details)

%plot 6 first spikes
spikeplot=[0 spiketrain(1:6)];
for i=1:length(spikeplot)-1
    hold on
plot(spikeplot(i+1),int(spikeplot(i+1)),'ko','MarkerFaceColor','k','MarkerSize',5,'LineWidth',2.5)
points=[spikeplot(i)+tau:0.0001:spikeplot(i+1)];
plot(points,int(points),'black','Linewidth',2);
end
shg

%Compute the log-likelihood analytically


loglik=log(1-cdftrans(T-spiketrain(end)));

for k=1:length(spiketrain)-1
    loglik=loglik+log(pdftrans(spiketrain(k+1)-spiketrain(k)));
end
loglik=loglik+log(pdftrans(spiketrain(1)));




disp('Spike Train Simulation Done')
disp('The plot contains the CIF during the first 6 spikes') 


%Number of evaluations
M=5000;

disp(['The true, analytically computed log likelihood, is loglik=' num2str(loglik)])
disp(['Now  approximations will be computed. By default, M=' num2str(M) ])
disp('log-likelihood evaluations')
disp('Press any key to compute the approximations')
pause

%minimum number of evaluations per integral
mineval=2;



t=cputime;
disp(['Now obtaining the Discrete Riemann 1 approximation'])
spikebins=BinSpikes(spiketrain,T,M);
points=linspace(0,T,M+1);
lags=tLastSpike(points(1:end-1),spiketrain);

loglikDR1=LoglikelihoodRiemannDerivatives(T,spiketrain,spikebins,lags,int,@(x)log(int(x)),tau,tau2,M,1);
disp(['DR1 approximation done. It took ' num2str(cputime-t) ' seconds to compute'])


t=cputime;
disp(['Now obtaining the Discrete Riemann 2 approximation'])
spikebins=BinSpikes(spiketrain,T,M);
points=linspace(0,T,M+1);
lags=tLastSpike(points(1:end-1),spiketrain);

loglikDR2=LoglikelihoodRiemannDerivatives(T,spiketrain,spikebins,lags,int,@(x)log(int(x)),tau,tau2,M,2);
disp(['DR2 approximation done. It took ' num2str(cputime-t) ' seconds to compute'])


t=cputime;
disp(['Now obtaining the Continuous Trapezoidal approximation'])

[node weight neval]=NodesWeightsCT(spiketrain,T,tau,M,mineval);
lagsCT=tLastSpike(node,spiketrain);
loglikCT=LoglikelihoodGLDerivatives(T,spiketrain,lagsCT,int,@(x)log(int(x)),tau,tau2,weight,node);
disp(['Gontinuous Trapezoid approximation done. It took ' num2str(cputime-t) ' seconds to compute'])


t=cputime;
disp(['Now obtaining the Gauss-Lobatto approximation'])

[node weight neval]=NodesWeightsGL(spiketrain,T,tau,M,mineval);
lagsGL=tLastSpike(node,spiketrain);
loglikGL=LoglikelihoodGLDerivatives(T,spiketrain,lagsGL,int,@(x)log(int(x)),tau,tau2,weight,node);
disp(['Gauss-Lobatto approximation done. It took ' num2str(cputime-t) ' seconds to compute'])
disp(['Summary:'])
disp(['True value of log likelihood    = ' num2str(loglik) ])
disp(['DR1 approximation is  loglikDR1 = ' num2str(loglikDR1)])
disp(['DR2 approximation is  loglikDR2 = ' num2str(loglikDR2)])
disp(['CT  approximation is  loglikCT  = ' num2str(loglikCT)])
disp(['GL  approximation is  loglikGL  = ' num2str(loglikGL)])
close(gcf)