function y=RenewalIntensity(t,f,spiketrain)
%Given an spike train driven by a renewal process
%with ISI intensity $f$
%creates the corresponding intensity function, to
%be evaluated at time t
for i=1:length(t)
t2=find(t(i)>spiketrain);
if(isempty(t2))

y(i)=f(t(i));
else

y(i)=f(t(i)-spiketrain(t2(end)));
end
end