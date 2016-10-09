function count=BinSpikes(spiketrain,T,M)

%Given the Spike Train spiketrain, recording time T and
%number of bins M, outputs the 0-1 time series count
%count(i)=1 if there is a spike, 0 if not

points=linspace(0,T,M+1);

for i=1:length(points)-1
    
    aux1=length(find(points(i+1)<=spiketrain));
    aux2=length(find(points(i)<=spiketrain));
    count(i)=aux2-aux1;
end
    