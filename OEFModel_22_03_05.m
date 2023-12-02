% This file contains the models of the group assignment "....". It consists of two scenarios,
% in which the cost of hydrogen production is calculating by assesing the CAPEX/X. 
%% Scenarios
% - Scenario On: Electricity is produced at an OWP, transported via HVDC to shore 
%   and converted into hydrogen via onshore electrolysis
% - Scenario Off: Hydrogen is production centralized on an offshore electrolysis platform for and then transported via pipeline to
%   shore
%% Model Structure
% I. Energy losses and demand
%   I.1 Scenario-independent calculations
%   I.2 Onshore Case  
%   I.3 Offshore Case
% II. Energy production
%   II.1 Weibull Distribution, Power Curve of OWP
%   II.2 Net power for electrolyzer, curltailment and partial-load
% III. Cost assessment
% IV. LCOH
% V. Analysis of model modifications
% VI. Plotting of relevant graphs
%% I. Energy losses and demand
%% I.1 Scenario-independent calculations
P_elec_A = zeros;
P_elec_B = zeros;
P_nom_elec = 100; % 100 MW rated power of the electrolyzer
eta_elec = 0.655; % conversion ratio of electricity into hydrogen
h2_endens = 120; % MJ/kg gravimetric energy density of hydrogen
h2_rho = 0.08988;% g/l density of hydrogen

E_h2_dmax = eta_elec*P_nom_elec*24; % MWh daily maximum of hydrogen production 
m_h2_dmax = (E_h2_dmax*3600) / h2_endens;%kg daily maximum of hydrogen production
m_h2_ymax = m_h2_dmax * 365;% kg yearly maximum of hydrogen production
vol_h2o_dmax = 10*m_h2_dmax/1000; % daily water volume demand, assuming 10kg water for 1 kg hydrogen, and 1000kg/m3 water density

dsal_el = 3;%kWh per m^3 per day electricity consumption of desalination unit
E_dsal = dsal_el*vol_h2o_dmax; %3-4kWh/m^3 elecricity demand for the desalination unit
P_dsal = 0.001*dsal_el*vol_h2o_dmax/24; %MW power capacity of desalination unit

comp_el = 12.47; %MJ/kg_h2 electricity consumption for compressor unit
eta_array = 0.9945; %efficiency of inter array cabling collecting AC from the turbines

%% I.2 Onshore Case  
eta_hvdc = 0.99*0.99*0.999965^140; % electricity losses in the substations and HVDC transmission to shore 
E_comp_A = comp_el*m_h2_dmax;%MJ electricity demand for the compressor
P_comp_A = E_comp_A/(24*60*60); %MW power capacity for the compressor
%% I.3 Offshore Case
eta_inv = 0.97; % inverter efficiency
eta_pipeline = 0.97; % pipeline efficiency
E_comp_B = comp_el*m_h2_dmax;%MJ electricity demand for the compressor
P_comp_B = E_comp_B/(24*60*60); %MW power capacity for the compressor
%% II Energy production
T = 15; % project duration, # years
totalh = 8760; % hours of one year
v = 1:25; % array for velocities of Weibull distribution
v = v.';
v_in = 3.5; %cut in speed
v_r = 14; %rated speed
v_out = 25; %cut out speed
d_tw = 193;% rotor diameter
P_wt = 10; % rated power of 1 turbine
N_wt = 12; % number of turbines
P_owp = P_wt * N_wt; % OWP capacity

P_gross = zeros(length(v),1);
P_net_A = zeros(length(v),1);
P_net_B = zeros(length(v),1);
P_curl_A = zeros(length(v),1);
P_curl_B = zeros(length(v),1);
P_pot_A = zeros(length(v),1);
P_pot_B = zeros(length(v),1);
E_elprod = zeros(length(v),1);
E_elec_A = zeros(length(v),1);
E_elec_B = zeros(length(v),1);
E_curl_A = zeros(length(v),1);
E_curl_B = zeros(length(v),1);
E_pot_A = zeros(length(v),1);
E_pot_B = zeros(length(v),1);



%% II.1 Weibull Distribution, Power Curve of OWP
h = zeros(length(v),1); % frequency according to Weibull
for i = 1:25
k = 2.3; % Weibull curve parameter
A = 10.8; % Weibull parameter for avg wind speed
h(i) = wblpdf(v(i), A, k);

% Power curve of wind park
% Taken from data sheet of a 10MW offshore turbine: 
% source: https://www.thewindpower.net/turbine_en_1662_siemens-gamesa_sg-10.0-193-dd.php
P_gross = N_wt*[0, 0, 0, 0.209, 0.801, 1.713, 2.92, 4.333, 6.201, 8.092, 9.299, 9.962, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]';

% Alternative code for approximating power curve from a cubic function
% We opted for the other method because its data provided by a manufacturer
% of a 10 MW turbine
% a= P_wt/(v_r^3-v_in^3); % parameter for cubic power curve
% b=-(P_wt*v_in^3)/(v_r^3-v_in^3); % parameter for cubic power curve
% if v(i) >= v_in && v(i) < v_r % partial load phase of power curve
%     P_gross(i) = N_wt*(a*v(i)^3+b);
% elseif v(i) >= v_r && v(i) < v_out % nominal load phase of power curve
%     P_gross(i) = P_owp;
% else 
%     P_gross(i)=0;  
% end

%% II.2 Net power for electrolyzer, curltailment and partial-load
P_net_A(i) = max(0, P_gross(i)*eta_array*eta_inv*eta_hvdc - P_dsal - P_comp_A);
P_net_B(i) = max(0, P_gross(i)*eta_array*eta_inv - P_dsal - P_comp_B);

if P_net_A(i) >= P_nom_elec % max load, curtailment
    P_elec_A(i,1) = P_nom_elec;
    P_curl_A(i) = P_net_A(i) - P_nom_elec;
else % partial load
    P_elec_A(i,1) = max(0,P_net_A(i));
    P_curl_A(i) = 0;
    P_pot_A(i) = P_nom_elec - P_net_A(i);
end 
if P_net_B(i) >= P_nom_elec
    P_elec_B(i,1) = P_nom_elec;
    P_curl_B(i) = P_net_B(i) - P_nom_elec;
else
    P_elec_B(i,1) = max(0,P_net_B(i));
    P_curl_B(i) = 0;
    P_pot_B(i) = 100 - P_net_B(i);
end
E_elec_A(i) = h(i)*P_elec_A(i)*totalh;
E_elec_B(i) = h(i)*P_elec_B(i)*totalh;
E_curl_A(i) = h(i)*P_curl_A(i)*totalh;
E_curl_B(i) = h(i)*P_curl_B(i)*totalh;
E_pot_A(i) = h(i)*P_pot_A(i)*totalh;
E_pot_B(i) = h(i)*P_pot_B(i)*totalh;
E_elprod(i) = h(i)*P_gross(i)*totalh;
end

P_owp_avg = sum(h.*P_gross); %MW
E_elprod_year = sum(E_elprod); %MWh

E_elec_A_year = sum(E_elec_A); %Mwh
E_elec_B_year = sum(E_elec_B); %MWh
E_curl_A_year = sum(E_curl_A);
E_curl_B_year = sum(E_curl_B);
E_pot_A_year = sum(E_pot_A);
E_pot_B_year = sum(E_pot_B);
m_h2_year_A = E_elec_A_year*3600/h2_endens;
m_h2_year_B = E_elec_B_year*3600/h2_endens;

%Capacity factors
cap_fac = E_elprod_year/(totalh*P_owp); % capacity factor
cap_fac_elec_A = E_elec_A_year/(totalh*P_nom_elec);
cap_fac_elec_B = E_elec_B_year/(totalh*P_nom_elec);
hours_elec_A = cap_fac_elec_A*totalh*T;
hours_elec_B = cap_fac_elec_B*totalh*T;

%% III. Cost Assessment 
d_shore = 140;
d_shore_analysis = 20:20:300;
%% Capex offshore wind park (OWP)
avg_cost_owp = 1000*4143*0.88; % average installed cost in USD/kW, converted into EUR, source: https://www.irena.org/publications/2021/Jun/Renewable-Power-Costs-in-2020
Cp_owp = avg_cost_owp * P_owp;
Op_owp = 0.0025 * Cp_owp;
%% CapEx electrolyzer and desalination
costpMW_elec = 650000; %€/MW, source: https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_renewable_fuels.pdf
Cp_elec = P_nom_elec * costpMW_elec;
Cp_desal = 0.88*1600*10*vol_h2o_dmax; % source: http://dx.doi.org/10.1016/j.desal.2012.10.015
%% CapEx Onshore Scenario (A)
Cp_sub_A = 193.24*N_wt*P_wt*1000; %€, cost for DC/DC substation, not including platform cost % http://dx.doi.org/10.3390/en12163191 
Cp_plat_A = 12000000; %€, estimated cost for platform of offshore substation, https://www.tennet.eu/fileadmin/user_upload/Company/Publications/Technical_Publications/Dutch/P2H_IJmuiden_Ver_-_Final_Report_-_Public.pdf
cost_hvdc = 2*760284; %€/km median cost for two HVDC subsea cable https://documents.acer.europa.eu/Official_documents/Publications/UIC_Electricity_History/UIC%20report%20%20-%20Electricity%20Infrastructure%20corrected.pdf
Cp_hvdc = cost_hvdc * d_shore;
Cp_comp_A = 2802000*E_comp_A/(24*60*60);%€ compressor cost calculated with 2802€ per kW source: https://projecten.topsectorenergie.nl/storage/app/uploads/public/5d0/263/410/5d026341016a2991247120.pdf
%% OpEx Onshore Scenario (A)
Op_sub_A = 0.002 * Cp_sub_A;
Op_plat_A = 0.002 * Cp_plat_A;
Op_hvdc = 0.0015 * Cp_hvdc;
Op_desal_A = 0.002 * Cp_desal;
Op_elec_A = 0.002 * Cp_elec;
Op_comp_A = 0.002 * Cp_comp_A;
%% CapEx Offshore Scneraio (B)
Cp_sub_B = 0.4*Cp_sub_A;
Cp_plat_B = 17000000;%€ estimate for a three level 30m x 35m platform based on existing oil&gas platform,https://www.tennet.eu/fileadmin/user_upload/Company/Publications/Technical_Publications/Dutch/P2H_IJmuiden_Ver_-_Final_Report_-_Public.pdf
cost_pipe = 1000000;%€/km source, p. 68: https://ore.catapult.org.uk/wp-content/uploads/2020/09/Solving-the-Integration-Challenge-ORE-Catapultr.pdf
Cp_pipe = cost_pipe*d_shore;
Cp_comp_B = 2802000*E_comp_A/(24*60*60);%compressor cost calculated with 2802€ per kW source: https://projecten.topsectorenergie.nl/storage/app/uploads/public/5d0/263/410/5d026341016a2991247120.pdf
%% OpEx Offshore Scneraio (B)
%The total OPEX are assumed to be about 2-3% of the total CAPEX.
%Over a project span of 15 years, the yearly OPEX accounts for
%0,2% to 0,3% of the CapEx, source: 10.1016/j.rset.2021.100005
Op_plat_B = 0.002 * Cp_plat_B;
Op_sub_B = 0.002 * Cp_sub_B;
Op_desal_B = 0.003 * Cp_desal;
Op_elec_B = 0.003 * Cp_elec;
Op_pipe = 0.002 * Cp_pipe;
Op_comp_B = 0.002 * Cp_comp_B;
%% Total CapEx and OpEx
df = 0.05; % discount factor
pvf = (((1+df)^T)-1)/(df*(1+df)^T); %pension present value factor

CapEx_A_vec = [Cp_owp Cp_plat_A Cp_elec Cp_desal Cp_sub_A Cp_hvdc Cp_comp_A];
CapEx_B_vec = [Cp_owp Cp_plat_B Cp_elec Cp_desal Cp_sub_B Cp_pipe Cp_comp_B];
CapEx_A = sum(CapEx_A_vec);
CapEx_B = sum(CapEx_B_vec);

OpEx_A_vec = [Op_owp Op_plat_A Op_elec_A Op_desal_A Op_sub_A Op_hvdc Op_comp_A];
OpEx_B_vec = [Op_owp Op_plat_B Op_elec_B Op_desal_B Op_sub_B Op_pipe Op_comp_B];
OpEx_A = sum(OpEx_A_vec);
OpEx_B = sum(OpEx_B_vec);

%% IV. Calculation of LCOH
LCOH_A_kg = (CapEx_A + pvf*OpEx_A)/(m_h2_year_A*pvf);
LCOH_B_kg = (CapEx_B + pvf*OpEx_B)/(m_h2_year_B*pvf);

LCOH_A_kWh = (CapEx_A + pvf*OpEx_A)/(m_h2_year_A*pvf)/(h2_endens/3.6);
LCOH_B_kWh = (CapEx_B + pvf*OpEx_B)/(m_h2_year_B*pvf)/(h2_endens/3.6);

%% V. Analysis of parameter modification
%% V.1 Distance to shore
Cp_hvdc_shore = zeros;
Cp_pipe_shore = zeros;
Op_hvdc_shore = zeros;
Op_pipe_shore = zeros;
CapEx_A_shore = zeros;
CapEx_B_shore = zeros;
OpEx_A_shore = zeros;
OpEx_B_shore = zeros;
LCOH_A_kg_shore = zeros;
LCOH_B_kg_shore = zeros;
LCOH_A_kWh_shore = zeros;
LCOH_B_kWh_shore = zeros;

for i = 1:length(d_shore_analysis)
    Cp_hvdc_shore(i) = cost_hvdc * d_shore_analysis(i);
    Cp_pipe_shore(i) = cost_pipe*d_shore_analysis(i);
    Op_hvdc_shore(i) = 0.002 * Cp_hvdc_shore(i);
    Op_pipe_shore(i) = 0.002 * Cp_pipe_shore(i);
    CapEx_A_shore(i) = Cp_owp + Cp_sub_A + Cp_plat_A + Cp_hvdc_shore(i) + Cp_desal + Cp_elec;
    CapEx_B_shore(i) = Cp_owp + Cp_plat_B + Cp_desal + Cp_elec + Cp_elec + Cp_pipe_shore(i);

    OpEx_A_shore(i) = Op_owp + Op_sub_A + Op_plat_A + Op_hvdc_shore(i) + Op_desal_A + Op_elec_A;
    OpEx_B_shore(i) = Op_owp + Op_plat_B + Op_desal_B + Op_elec_B + Op_elec_B + Op_pipe_shore(i);

    LCOH_A_kg_shore(i) = (CapEx_A_shore(i) + pvf*OpEx_A_shore(i))/(m_h2_year_A*pvf);
    LCOH_B_kg_shore(i) = (CapEx_B_shore(i) + pvf*OpEx_B_shore(i))/(m_h2_year_B*pvf);

    LCOH_A_kWh_shore(i) = (CapEx_A_shore(i) + pvf*OpEx_A_shore(i))/(m_h2_year_A*pvf)/(h2_endens/3.6);
    LCOH_B_kWh_shore(i) = (CapEx_B_shore(i) + pvf*OpEx_B_shore(i))/(m_h2_year_B*pvf)/(h2_endens/3.6);
end

%% V.2 Analysis hydrogen pipeline cost per km
    cost_pipe_analysis = 500000:200000:1700000;
    Cp_pipe_analysis = zeros;
    Op_pipe_analysis = zeros;
    CapEx_B_pipe = zeros;
    OpEx_B_pipe = zeros;
    LCOH_B_kg_pipe = zeros;
 for i = 1:length(cost_pipe_analysis)
    Cp_pipe_analysis(i) = cost_pipe_analysis(i)*d_shore;
    Op_pipe_analysis(i) = 0.002 * Cp_pipe_analysis(i);
    CapEx_B_pipe(i) = Cp_owp + Cp_plat_B + Cp_desal + Cp_elec + Cp_elec + Cp_pipe_analysis(i);
    OpEx_B_pipe(i) = Op_owp + Op_plat_B + Op_desal_B + Op_elec_B + Op_elec_B + Op_pipe_analysis(i);
    LCOH_B_kg_pipe(i) = (CapEx_B_pipe(i) + pvf*OpEx_B_pipe(i))/(m_h2_year_B*pvf);
 end
clear i;
%% VI Plotting of relevant graphs
figure
hold on
box on
barh([CapEx_A_vec; CapEx_B_vec], 'stacked');
title('CAPEX of both Scenarios');
legend('OWP','Platform','Electrolyzer','Desalination','Substation', 'Transmission','Compressor');
xticks([0 1e8 2e8 3e8 4e8 5e8 6e8 7e8 8e8])
xlabel('CAPEX [M€]');
xticklabels({'0','100','200','300','400','500','600','700','800'})
yticks([1 2])
yticklabels({'A','B'})
hold off
figure
hold on
box on
barh([T*OpEx_A_vec; T*OpEx_B_vec], 'stacked');
title('OPEX of both Scenarios');
legend('OWP','Platform','Electrolyzer','Desalination','Substation', 'Transmission','Compressor');
xlabel('OPEX [M€]');
xticks([0,5000000,10000000,15000000,20000000,25000000,30000000])
xticklabels({'0','5','10','15','20','25','30'})
yticks([1 2])
yticklabels({'A','B'})
hold off
figure
hold on
box on
plot(d_shore_analysis,LCOH_A_kg_shore);
plot(d_shore_analysis,LCOH_B_kg_shore);
title('LCOH for different distances to shore');
xlabel('Shore distance [km]');
ylabel('LCOH [€/kg]');
legend("Onshore (A)","Offshore (B)");
hold off
figure 
hold on
box on
LCOH_pipe_vector(1,1:7) = LCOH_A_kg;
LCOH_pipe_vector(2,1:7) = LCOH_B_kg_pipe;
bar(cost_pipe_analysis, LCOH_pipe_vector);
title('LCOH comparison for varying pipeline cost')
legend('Onshore (A)','Offshore (B)');
ylabel('LCOH [€/kg]');
xlabel('Pipeline cost [M€/km]');
xticks([0.5e6 0.7e6 0.9e6 1.1e6 1.3e6 1.5e6 1.7e6]);
xticklabels({'0.5','0.7','0.9','1.1','1.3','1.5','1.7'});
hold off
% clearvars -except LCOH_A_kWh LCOH_B_kWh LCOH_A_kg LCOH_A_kg LCOH_B_kg