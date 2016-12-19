
Temperature_in_Kelvin = 300;
[t_sc,B] = inv_cm_unit_sys_standard(Temperature_in_Kelvin); %B = 1/kB T

mode_inc = 3;
switch mode_inc
    case 0 %all
lg = zeros(48,1);    
    case 1 %based on if BOmodes_tmp(:,3)>500 & BOmodes_tmp(:,1)<15 then exclude
lg = [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,0,1,0,1,1,1,0,0,1,0,0,1,1,1,1,1 ...
        ,0,1,1,1,1,1,1,0,1,1,1,1,1,1]'; %value of 1 removes from HEOM
    case 2 %BOmodes_tmp(:,3)>600 | BOmodes_tmp(:,1)<20
lg = [1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 ...
      ,1,1,1,1,1,1,1,1,1,1,1,1,1,1]';   %takes only 2 modes     
    case 3 %to take all out   
lg = ones(48,1) ; 
end
lg2 = true(48,1); %local vector to select particular modes

max_tier = 3;
Kappa = 1; %without high freq modes at room temp Kappa = 1 is fine
specden;  %call to get all the info including spec den and HEOM and Hamiltonian
plot_stuff = true; fit_rates = false;
[M_e,E1] = eig(H_1); %eigenstate basis

%% Generate full prop op for the dynamics
V_coup_tot = zeros(10,10,N); V_coup_sys = zeros(N,N,N);
for j = 1:N %project each term into the exciton basis
    tmp = zeros(N); tmp(j,j) = 1;    V_1 = M_e'*tmp*M_e;
V_coup_sys(:,:,j) = V_1; %zeros are g, alpha, beta, I
end                

prop_op_sys =H_prop_gen3(E1,E1,V_coup_sys,V_coup_sys,...
                QQ_topass,const_factor,coup_com_save,coup_acom_save);            
            
%prop_op_full(1:100,1:100) = prop_op_full(1:100,1:100)+L_total; %HEOM + rate
%this is one mother of a matrix!

%% Calculate system only dynamics

rnj = 1:N^2; pop_lg = false(N^2,1); pop_lg(1:N+1:N^2)=true;

[rho_ss,test] = eigs(prop_op_sys,1,'SM');
if abs(test) > eps(length(rho_ss)); warning('steady state not converged'); end
rho_ss = real(rho_ss / sum(rho_ss(pop_lg)));
rho_ss2 = reshape(rho_ss(rnj),N,N); rho_ss2 = (rho_ss2+rho_ss2')/2;
[bath_basis,pops] = eig(rho_ss2); %basis bath projects into
energies_sysonly = real(diag(bath_basis'*E1*bath_basis));
de_sys = energies_sysonly - max(energies_sysonly); 

%ff = @(E) exp(-B*E)./sum(exp(-B*E))-diag(pops); %gibbs eq
ff2 = @(E) [exp(-B*E);1]./(1+sum(exp(-B*E)))-diag(pops); %sets lowest to zero

%options.MaxFunEvals = 1000*N;
%[EE,FVAL] =fsolve(@(E) ff(E),energies_sysonly);
[EE2,FVAL2] =fsolve(@(E) ff2(E),de_sys(1:end-1)); %works better
tmp = diag(exp(-B*[EE2;0]));tmp = tmp/trace(tmp);
test = diag(tmp)'-diag(pops)'; %should be small
%EE2 are the energy differences (assuming a Gibbs state is accurate) of the
%states the bath projects into 

%difference between gibbs state
rho_gibbs = diag(exp(-B*diag(E1))); 
rho_gibbs = reshape(rho_gibbs/trace(rho_gibbs),N^2,1);
diff_mat = rho_gibbs - rho_ss(rnj);
trace_dist = sum(abs(eig(reshape(diff_mat,N,N))));

propeq_sys = @(t,v)  prop_op_sys*v;

rho_vec = zeros(N^2*size(nn,1),1);  rho_vec(1)=1;
t_range = linspace(0,1,100); 

[tout,rho_out] = ode45(propeq_sys,t_range,rho_vec);
%%
if plot_stuff
figure; plot1 = plot(tout,rho_out(:,pop_lg)); 
for k =1:sum(pop_lg)
set(plot1(k),'DisplayName',strcat('exciton ',num2str(k)));
end
xlabel('time (ps)'); ylabel('exciton population')
%%
tmp = reshape(rho_out(:,1:length(pop_lg)).' ,[N,N,length(tout)]);
coherence_set = zeros(sum(~pop_lg)/2,length(tout)); cnt = 0;
for k1 = 1:N; for k2=k1+1:N;
cnt = cnt+1;
coherence_set(cnt,:) = tmp(k1,k2,:)./sqrt(tmp(k1,k1,:).*tmp(k2,k2,:));
end; end     
        
figure; plot2 = plot(tout,abs(coherence_set)); cnt=0;
for k1 = 1:N; for k2=k1+1:N;
cnt = cnt+1;
set(plot2(cnt),'DisplayName',strcat('\rho_{',num2str(k1),num2str(k2),'}'));
end; end     

xlabel('time (ps)'); ylabel('normalised coherence $|\rho_{ab}|/\sqrt{\rho_{aa}\rho_{bb}}$')

end
%% Fit rates
if fit_rates 
pop_t = rho_out(:,pop_lg).';
[R_ab,minval,exitflag] = fit_transfer_rates(pop_t,tout);

lg = logical(reshape(eye(N),[N^2,1]));
prop_mat = zeros(N); prop_mat(~lg) = R_ab; 
prop_mat(lg) = -sum(prop_mat,1); 
f = @(t,x) prop_mat*x;
[~,xs] = ode45(f,tout, pop_t(:,1)); 
figure; plot(tout,xs); hold on; plot(tout,pop_t.','--')
end
