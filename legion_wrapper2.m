%legion_wrapper edit to check convergence

max_tier = 10; 
Kappa = 0;
Temperature_in_Kelvin = 300;
reorg_energy_rnj = [35,70,100,150];  %range of reorganisation energies
t_range = linspace(0,0.1,200); %range of time to prop for
save_tier_cut = 1;
%% Don't edit the rest

lg = ones(48,1) ; 
lg2 = true(48,1); %local vector to select particular modes

%%
for some_lp = 1:length(reorg_energy_rnj )
    new_rorg = reorg_energy_rnj(some_lp);
alternative_reorg = [new_rorg,40];
specden;  %call to get all the info including spec den and HEOM and Hamiltonian
[M_e,E1] = eig(H_1); %eigenstate basis       
    
% Generate full prop op for the dynamics
V_coup_tot = zeros(10,10,N); V_coup_sys = zeros(N,N,N);
for j = 1:N %project each term into the exciton basis
    tmp = zeros(N); tmp(j,j) = 1;    V_1 = M_e'*tmp*M_e;
V_coup_sys(:,:,j) = V_1; %zeros are g, alpha, beta, I
end                
prop_op_sys =H_prop_gen3(E1,E1,V_coup_sys,V_coup_sys,...
                QQ_topass,const_factor,coup_com_save,coup_acom_save);            

propeq_sys = @(t,v)  prop_op_sys*v;

rho_vec = zeros(N^2*size(nn,1),1);  rho_vec(1)=1;

[tout,rho_out] = ode45(propeq_sys,t_range,rho_vec);
flename = strcat('exciton_dynamics_reorg',num2str(new_rorg),'MT',...
                num2str(max_tier),'Kap',num2str(Kappa),'.mat');

%save only part of rhoout, to whatever tier you decide to cut it off
cutoff = sum(nn,2)<=save_tier_cut;  %for what is saved
rnj = 1:N^2*sum(cutoff); rho_out = rho_out(:,rnj);

save(flename,'tout','rho_out','max_tier','Kappa','new_rorg')

end