% Input parameters for Cees Dekker, PNAS, SiNx ion transport 
% measurement data 2006. % KCl ionic solution

% Figure 2(a) conc = 100 mM, surface charge density = - 0.5 mC/m2
% As expected Bulk dominates at this conc, sigma_expt = -20 mC/m2, but our sigma not equal to
% paper fit results. A direct evidence of sigma is needed to complete the
% theory


clear;
clc;

% Simulation parameters for single potassium ion transport
D = 1.96e-9 ; % Diffusivity of Potassium
K_b = 1.38064e-23; % in kg m^2/s^2K Boltzman Constant
T = 300; % Temperature in K
z = 1; % Valence of Potassium ion
lambda = (K_b*T)/D ; % Mass rate in (kg/s)
omega = D/(K_b*T) ; % Mobility in (s/kg)
epsilon_0 = 8.854e-12; % Dielectric permittivity of free space in F/m
epsilon_r = 80; % Dielectric permittivity of water (no units)
e_charge = 1.60217e-19; % electronic Charge

% Geometry of the nanoporous membrane

r_nano = 5e-9; % Radius of nanopore in nm, diameter - 10 nm for calculation
L_nano = 34e-9; % Length of nanopore in nm, length = 34 nm for calculation

A_nano = (pi*(r_nano)^2); % cross sectional-area of nanopore

dS_nano = (2*pi*r_nano*L_nano); % Surface area of nanopore
dV_nano = A_nano*L_nano; % Volume of the nanopore

sigma_nano = -0.5e-3; % surface charge of nanopore in C/m^2, negative surface charge density




%%%%%%%%%%%%%%%%%%%% concentration of the ionic solution %%%%%%%%%%%%%%%

c = 100; % concentration of KCl = 100 mM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NA = 6.023e23; % Avagadro number

% Total mass of the Charged particle
mass_Potassium = 39.10e-3; % Kg/mol
mass_single_atom = ((mass_Potassium)/NA); % mass of one charged ion


% Simulation time step
t_new = (1e-15:1e-15:1e-12); % time considered in fs to ps range with dt=1fs
N_size = size(t_new,2);% number of time steps

% Applied Input voltage
Voltage = (-34e-3:1e-3:34e-3) ;
size_V = size(Voltage,2);

% Derived parameters
p = 4*pi*epsilon_0*epsilon_r; % constant
x_b1 = 0;
x_b2 = L_nano;
L_total = L_nano;


% Initial guess of each parameter
q_new = ones(1,N_size)*(e_charge*z); % Initial total counter ionic charge (GUESS 1)
x_variable(1) = 1e-9; % initial position of the ion
v_new(1) = 0.2; % initial velocity = ~ intial displacement (= 1nm)/diffusion time (= 1ns) in m/s
n_total = (c*NA)*(dV_nano); % total number of Potassium ions
q_s_nano = sigma_nano*(2*pi*r_nano*L_total); % charge of nanopore
current(1) =(v_new(1)*q_new(1)*n_total)/L_total;

%To calculate individual forces at time t = 1fs
F1(1) = (abs(q_new(1))*Voltage(1))/((x_b2-x_b1));
F2(1) = abs(q_new(1))*abs(q_s_nano)/(4*pi*epsilon_0*epsilon_r*((L_nano-x_variable(1))^2+(r_nano)^2));
F3(1) = (abs(q_s_nano)*Voltage(1))/((x_b2-x_b1));
Fdissipative(1) = (lambda*v_new(1));

% acceleration at time step 1

acc_new(1) = F1(1)/mass_single_atom + F2(1)/mass_single_atom + F3(1)/mass_single_atom - Fdissipative(1)/mass_single_atom ; 

for j = 1:size_V
    
for i = 2:N_size % From time step 2 to N_size; number of time steps

% position of the particle calculated using Velocity-Verlet Algorithm

x_variable(i) = x_variable(i-1) + (v_new(i-1)*(t_new(i)-t_new(i-1))) + (0.5*(acc_new(i-1)*(t_new(i)-t_new(i-1))^2));

%To calculate individual forces at time = t 
F1(i) = (abs(q_new(i))*Voltage(j))/((x_b2-x_b1));
F2(i) = abs(q_new(i))*abs(q_s_nano)/(4*pi*epsilon_0*epsilon_r*((L_nano-x_variable(i))^2+(r_nano)^2));
F3(i) = (abs(q_s_nano)*Voltage(j))/((x_b2-x_b1));
Fdissipative(i) = (lambda*v_new(i-1));

% acceleration of the particle
acc_new(i) = F1(i)/mass_single_atom + F2(i)/mass_single_atom + F3(i)/mass_single_atom - Fdissipative(i)/mass_single_atom ; 
%velocity of particle
v_new(i) = v_new(i-1) + (0.5*(acc_new(i-1)+acc_new(i))*(t_new(i)-t_new(i-1)));
% Current in the nanopore
current(i) = (v_new(i)*q_new(i)*n_total)/L_total;
   
  
 if (x_variable(i) > L_nano)
     size_plot = i;
  %  break;
 end 
 
  if (x_variable(i) < 0)
     size_plot = i;
  %   break;
  end
  
  size_plot = i;
end

CURRENT(j,1) = current(size_plot)
conductance(j,1) = CURRENT(j,1)/Voltage(j);
end

% Tabulated values
Voltage_Dekker_PNAS = Voltage'*1e3; %(mV)
CURRENT_Tabulate = CURRENT*1e12;% (pA)

hold on; plot(Voltage*1e3, CURRENT*1e12, '-m');
xlabel('Voltage (mV)', 'fontsize', 24);
ylabel('I (pA)', 'fontsize', 24);

