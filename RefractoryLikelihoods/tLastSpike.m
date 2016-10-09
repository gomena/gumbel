function y=tLastSpike(t,spiketrain)

for i=1:length(t)
    ind=find(t(i)>spiketrain);
    if(length(ind)>0)
    y(i)=t(i)-spiketrain(ind(end));
    else
        y(i)=t(i);
    end
end