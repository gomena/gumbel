function x=NewtonBackTrack(f,df,d2f,tol,maxiter,x0)

%minimizes f with the Newton Method with backtracking line-search
%df and df2 are the gradient and hessian, respectively
%x0 is the initial point
%repeat until the norm of the gradient is smaller than tol
%or until exceding a maximum number of iterations, maxiter
cont=1;
x=x0;

c=0.1;
alpha0=1;
tau=0.5;
grad=df(x);
cont=1;
while(norm(grad)>tol&&cont<=maxiter)
    

    alpha=alpha0;
    dir=-d2f(x)\grad;
    t=dir'*grad*c;
    while(f(x)>f(x+alpha*dir)-alpha*t)
        alpha=alpha*tau;
    end
    
    x=x+alpha*dir;
    grad=df(x);
    
  disp(['Newton Iteration ' num2str(cont) ' finished. Norm of the gradient = ' num2str(norm(grad))])
  cont=cont+1;
end
if(cont>maxiter)
    disp(['Maximum number of iterations exceeded'])
end
        
    
    