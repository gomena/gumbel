function [nodes weights]=CreateNodes(n,type)

%create the first n+1 nodes and weights vectors using the
%Gauss-Lobato-Legendre rule (type=1) or the Trapezoidal rule (type=2) 
%the vector nodes{i} contains the nodes of the interval 
%using i+1 points
%the vector weights{i} contains nodes of the interval 
%using i+1 points
if(type==1)
for i=1:n
    
    [x,w,P]=lglnodes(i);
    nodes{i}=x;
    weights{i}=w;
end
else
    for i=1:n
        nodes{i}=linspace(-1,1,i+1)';
        weights{i}=2/i*ones(i+1,1);
        weights{i}([1 end])=1/i;
    end
end
    