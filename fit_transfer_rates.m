function [R_ab,minval,exitflag] = fit_transfer_rates(Pop_t,trange,p0)
%rate fitting code:  Fits set of transfer rates to a time dependent data
%set to obtain rates of transfer from state a into state b
% Pop_t are time dep population rates
% p0 is an initial guess (based on Markovian rates / Modified Redfield etc

TL = length(trange);
N = size(Pop_t,1); 
if N == TL && size(Pop_t,2) ~= TL; N = size(Pop_t,2); Pop_t=Pop_t.';  end
if nargin < 3
    p0 = rand(N*(N-1),1); 
end

init_state = Pop_t(:,1);
options.MaxFunEvals = 1000*length(p0);

 [R_ab , minval, exitflag] = fminsearch(@(p) myobjective(trange,Pop_t,p,init_state), p0,options);
 
end

function prop_mat = prop_mat_gen(p,N)

    lg = logical(reshape(eye(N),[N^2,1]));
    prop_mat = zeros(N); 

    prop_mat(~lg) = p; 
    prop_mat(lg) = -sum(prop_mat,1);        

end


function SSE = myobjective(td,xd, p, x0)

        prop_mat = prop_mat_gen(p,length(x0));
         f = @(t,x) prop_mat*x;  %function to call ode45
         [~,xs] = ode45(f, td, x0);  %xs is the state variable predicted by the model
         if size(xs)~=size(xd); xd = xd.'; end
         err = xd - xs;
         SSE = sum(sum(err.^2));   %sum squared-error. 
         SSE = SSE + sum(abs(p(p<0))); %negative weights are baddd
end