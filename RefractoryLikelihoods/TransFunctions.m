function [pdftau,cdftau,inttau]=TransFunctions(pdf,cdf,int,tau)

%Given the pdf,cdf and intensity functions of the ISI
%returns their tau-translated version (see paper)

pdftau=@(x)pdf(x-tau).*(x>=tau);
cdftau=@(x)cdf(x-tau).*(x>=tau);
inttau=@(x)int(x-tau).*(x>=tau);