function [x w neval]=NodesWeightsCT(varargin)

%this .m file does essentialy the same as NodesWeightsGL.m
%but it loads the nodes and weights using the Trapezoidal rule instead
spiketrain=varargin{1};
T=varargin{2};
tau=varargin{3};
M=varargin{4};
if(nargin>4)
mineval=varargin{5};
else
    mineval=1;
end
if(nargin>5)
    p=varargin{6};
else
    p=100;
end
%Computes the evaluation points(x) in the interval [0,T]%
%and the corresponding  weights(w) using the Gauss-Lobatto
%quadrature rule for each ISI. The number of evaluation points
%per ISI is given in the vector neval

%Arguments
%spiketrain=vector of spiking times
%T=Recording time
%tau=Absolute refractory period 
%M=number of evaluation points to be distributed in the ISI's
%mineval=Minimum of evaluation points to be assigned to each integral
%(by default mineval=1)
%p= Optional parameter. Finds pr, the p-percentile of the
%distribution of evaluation points among the ISI's
%then all the ISI's with more than pr evaluation points are
%reduced to pr, and the remaining evaluation points are
%redistributed randomly into the others.
%by default p=100, so this operation is not done.

load nodesweightsCT
nmax=size(weights,2);
%Maximum of trapezoid rule nodes to be used in each integral.
%nmax=300 by default. To try different values create a different
%set of nodes and weights and save it to the file nodesweightsCT.mat
%(to this end use the function CreateNodes.m)
%If the number of evaluations exceeds nmax then the interval is subdivided
%and the nodes and weights are computed for each of the sub-divided
%intervals



nint=length(spiketrain)+1;
flag=1;
if((mineval)*nint>M)
 
   
disp('Not enough evaluation points.' )
disp(['M was set as the minimum possible to' ...
'make all the required computations, M=' num2str(mineval*nint)] )
    flag=0;
   

end
%Now the number of evaluations per integral
%Distribute the remaining points proportionally to the length
%of each ISI


T2=T-(nint-1)*tau;    
neval=mineval*ones(1,nint);
M2=M-mineval*nint-1;

staux=[0 spiketrain T];
prop=floor(M2/T2*(staux(2:end)-staux(1:end-1)-tau));
prop(1)=floor(M2/T2*(spiketrain(1)));
neval=neval+prop;

rem=M-sum(neval);
if(rem>0)
add=mnrnd(rem,ones(nint,1)/nint);    

neval=neval+add;
end

%Redistribute points such that none of the ISI's receives
%more than the p percentile.

pr=floor(prctile(neval,p));
dif=sum(max(0,neval-pr));
neval(neval>pr)=pr;
if(dif>0)
  
add=mnrnd(dif,ones(sum(neval<=pr),1)/sum(neval<=pr));
neval(neval<=pr)=neval(neval<=pr)+add;
end
spiketrainext=[0 spiketrain T];


if(flag==0);
    neval=mineval*ones(nint,1);

end
neval(1)=max(neval(1),2);


%Now find nodes are weigths for each ISI 
x=[];w=[];
for i=1:nint
    
    if(neval(i)<=nmax)
       
        
        if(i==1)
%First ISI is different as there is no assumed refractoriness

x=nodes{neval(1)-1}'*(spiketrain(1)/2)+spiketrain(1)/2;
x(end)=spiketrain(1);
w=spiketrain(1)/2*weights{neval(1)-1}';
        else

x=[x nodes{neval(i)}(2:end)'*(spiketrainext(i+1)-spiketrainext(i)-tau)/2+(spiketrainext(i+1)+spiketrainext(i)+tau)/2];
x(end)=spiketrainext(i+1);

w=[w weights{neval(i)}(2:end)'*(spiketrainext(i+1)-spiketrainext(i)-tau)/2];

        end
    else
%If number of evaluations exceeds nmax then subdivide the interval in
%the minimum number, nint2, of subintervals with equal lenght such that 
%the number of evaluations in each sub-interval doesn't exceed nmax.

    nint2=floor((neval(i)-1)/nmax)+1;
    
    neval2=[floor(neval(i)/nint2)*ones(1,nint2)];
    if(sum(neval2)<neval(i))
        rem=neval(i)-sum(neval2);
        neval2(1:rem)=neval2(1:rem)+1;
    end
    
    
    if(i==1)
     %Again, the first ISI is treated slightly differently
        interval=linspace(spiketrainext(i),spiketrainext(i+1),nint2+1); 
     
    xaux=[]; waux=[];
   for k=1:nint2
   
xaux=[xaux nodes{neval2(k)-1}'*(interval(k+1)-interval(k))/2+(interval(k+1)+interval(k))/2]; 
waux=[waux weights{neval2(k)-1}'*(interval(k+1)-interval(k))/2];
   end
xaux(end)=spiketrainext(i+1);
xaux(end)=spiketrainext(i+1);
x=[x xaux];
w=[w waux];
    else
%When the interval is subdivided, the first subdivision (containing)
%the absolute refractory effects, is treated slighly different, to don't

    
  interval=linspace(spiketrainext(i)+tau,spiketrainext(i+1),nint2+1);  
    
    xaux=[];waux=[];
    xaux=[xaux nodes{neval2(1)}(2:end)'*(interval(2)-interval(1))/2+(interval(2)+interval(1))/2];
    waux=[waux weights{neval2(1)}(2:end)'*(interval(2)-interval(1))/2]; 
    
    for k=2:nint2
        
    xaux=[xaux nodes{neval2(k)-1}'*(interval(k+1)-interval(k))/2+(interval(k+1)+interval(k))/2];
    waux=[waux weights{neval2(k)-1}'*(interval(k+1)-interval(k))/2]; 
    
    end

xaux(end)=spiketrainext(i+1);
xaux(end)=spiketrainext(i+1);
x=[x xaux];
w=[w waux];
    end
    end
end