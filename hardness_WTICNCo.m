
%% CCoNTiW system in 3 phase region at 1000C. 
%% Phases are FCC_A!#1(Co reach binder),FCC_A!#2(Gamma), MC_SHP(Alpha)
%% Components are Va C Co N Ti W
%% ELEMENTS are 1-C 2-Co 3-N 4-Ti 5-W
%% Model is (Co,W,Ti)1(C,N,Ti)1
%% Duo to lack of data system is equilalized to CCoW with only two phase: 
%% MC_SHP, FCC_A1#1 binder
%% clear unnecessary variables
close all;
disp('!!!!!!!!!!!!<<hardness calc Start>>!!!!!!!!!!!!')
disp(work_path)
%% Open output files
fileID_1 = fopen([work_path '\HCC_HV.TXT'],'w');
fileID_2 = fopen([work_path '\K1C_Pa.TXT'],'w');
fileID_3 = fopen([work_path '\HCC_GPa.TXT'],'w');
fileID_4 = fopen([work_path '\K1C_MPa.TXT'],'w');
%% Dummy variables
ngp = NUMBER_OF_GRID_POINTS;
nel = NUMBER_OF_ELEMENTS;
ndim = DIMENSIONALITY;
dmsize = DOMAIN_SIZE;
mf = MOLE_FRACTIONS;
%% Calculate equilibrium in all points
%% creat a mole fraction condition string
%% Initialize TCMATLAB toolbox
%% Read a POLY3 file of the system
%% Set conditions & calculate equilibrium
%% Read phase names and number of phases
sel =[];
for i = 1 : nel
    X_of_el(i) = cellstr([ 'X(' ELEMENT_NEMAES{i} ')']);
    if i==1
        sel = strcat(sel,ELEMENT_NEMAES{i});
    else
        sel = strcat(sel, {' '}, ELEMENT_NEMAES{i});
    end
end
%% General conditions Thermodynamics
N = 1; % number of moles of the system
P = 101325; % Ambient pressure
T = 1273.15; % Freezed in T
%%
tc_init_root;
tc_read_poly3_file('wccotin_mcshp-fcca1-tcfe8.POLY3')
[nph_nw, ph_name_nw] = tc_list_phase;
[ncomp_nw, comp_name_nw] = tc_list_component;
tctoolbox('tc_set_condition', 'n', N);
tctoolbox('tc_set_condition', 'p', P);
tctoolbox('tc_set_condition', 't', T);
%% select the time step
tstp = nts;
%% Extract mole fractions of the timestep without baoundaries
%% this will be taken care of in the python postprocessor
if ndim == 2
    mf_tstp = mf((tstp-1)*nel*ngp(1)*ngp(2)+1: tstp*nel*ngp(1)*ngp(2));
end

%% General conditions Hardness
D_size = 1e-6; % grain size
k = 0.76; % for pure cobalt binder
R = 8.314; % Universal constant of Gases
TH = 298.15; % Hardness calculation T
%% One can make a matrix of inputs and solve it in a systematic way but in
%% this case the aproach is rather practical Va Co Ti W C N
H0Co = -107.5* log(TH)+ 789.48; 
H0W = -0.000000481* TH^3+ 0.00158*TH^2 + 736;
H0Ti = -3.44756259634154E-07*TH^3 + 0.00117493950541393*TH^2 - ...
    1.35042165191535*TH + 532.37586432309;
A_CoTi=56168;
A_CoW=4767;
A_TiW=239;
A_CoC=3354;
A_WC=0;
%A_TiC=A_WC;
A_TiC = 2000;
%A_CoN A_TiC A_TiN A_WN;   
Q_CoTi= 1024;
Q_CoW=280;
Q_TiW=298;
Q_CoC=2197;
Q_WC=294;
%Q_TiC = Q_WC ;  
Q_TiC = 200;
%Q_WN Q_CoN Q_TiC Q_TiN;
%% General conditions Toughness
    EB = 207*10^9; % Young's mpdulous of the Cobalt binder GPa
    EWC = 700*10^9; %Young's mpdulous of the WC GPa
    nuB = 0.3;% Poisson's ratio of the binder
    sigmaYB = 390*10^6;%Yield stress of the Co binder MPa 
    x0 = 250e-9;% initial void spacing 200-300 nm
    nuCC = 0.24;% poisson's ratio of cemented carbide
    GWC = 8;%fracture energy of WC 7-10 J/m^2
%     ECC = 5.6*10^11; % Pa from wiki 
%% loop over all nods - linear array
for i = 1 : size(mf_tstp,1)/nel
    f=1;
    %% Set concentration conditions
    for j = ((i-1)*nel+1): (i*nel-1)
        tc_check_error;
        tctoolbox('tc_set_condition',char(X_of_el(f)), mf_tstp(j));
        f = f+1;
    end
    %% Calcualte equilibrium
    tc_compute_equilibrium;
    %% Read the equilibrium % read VPV and Y fraction of Components of the Binder phase
    %(char(['vpv(' char(ph_name_nw{3}) ')']))
    vpvB = tc_get_value(char(['vpv(' char(ph_name_nw{1}) ')'])); % only Binder = fcc_a1#1
    vpvBinder(i) = vpvB;
    for comp_idx  = 1: ncomp_nw % only Binder = fcc_a1#1
        yf(comp_idx) = tc_get_value(char(['y(' char(ph_name_nw{1}) ',' char(comp_name_nw(comp_idx)) ')']));
    end
    YVa = yf(1);     
    YC  = yf(2) + yf(4); %YC == YC+YN
    YCo = yf(3);
    YTi = yf(5);
    YW  = yf(6);
    %% Hardness calculation Walbruhl's model
    H0 = H0Co*YCo*YVa + H0Ti*YTi*YVa+ H0W*YW*YVa; %(i)
    H_SSH = (1/9.807)*( ...
        A_CoTi * exp(Q_CoTi /R/TH) * (YCo* YTi)^2/3 * YVa + ...
        A_CoW  * exp(Q_CoW  /R/TH) * (YCo* YW )^2/3 * YVa + ...
        A_TiW  * exp(Q_TiW  /R/TH) * (YTi* YW )^2/3 * YVa + ...
        A_CoC  * exp(Q_CoC  /R/TH) * (YVa* YC )^2/3 * YCo + ...
        A_WC   * exp(Q_WC   /R/TH) * (YVa* YC )^2/3 * YW  + ...
        A_TiC  * exp(Q_TiC  /R/TH) * (YVa* YC )^2/3 * YTi ); %(i)
    HB  = H0 +  H_SSH; %(i)
    mfp = D_size*vpvB/ (1-vpvB); %(i)
    HCC(i) = ( (693+ 2680/sqrt( 2.1+ D_size))- HB) * exp(-((mfp*1e6)^0.67)/k) + HB; %(i) in Vickers / HV to GPa multiply by 0.009807
    %% Write Hardness to file
    fprintf(fileID_1, '%f \n', HCC(i)); %(i) in Vickers to GPA => HV multiply by 0.009807
    fprintf(fileID_3, '%f \n', HCC(i)*0.009807); %(i) in Vickers to GPA => HV multiply by 0.009807
    %% Fracture toughness calculation Linder WPMA 2018
    %GB = (sigmaYB*x0^(0.5)*(19.655+5.8508*log(mfp * 1e-6 /x0)))^2*(1-nuB^2)/EB;% fracture eenrgy of the binder
    GB = (1E6 *( 19.655+5.8508*log(mfp/x0) ))^2 * (1-nuB^2) / EB;% fracture eenrgy of the binder
    ECC = EB*((EB+(EWC-EB)*(1-vpvB).^(2/3))/(EB+(EWC-EB)*(1-vpvB).^(2/3)*(1-(1-vpvB).^(1/3)))); % Pa % EB is E_Co% %Paul relation Young modulus of cemented carbide
    K1C(i) = sqrt(  (ECC/(1-nuCC^2)) * (exp(-1.77*vpvB^0.78) * GWC + (1-exp(-1.77*vpvB^0.78))*GB)); %in Pa
    %% Write Thoughness to file
    fprintf(fileID_2, '%f \n', K1C(i)); %(i) in Pa
    fprintf(fileID_4, '%f \n', K1C(i)*1e-6); %(i) in MPa
end
%% Reshape linear arrays and surface 
HCC_2D = reshape( HCC, [ngp(1), ngp(2)] );
K1C_2D = reshape( K1C, [ngp(1), ngp(2)] );
vpv_2D = reshape( vpvBinder, [ngp(1), ngp(2)] );  
x = FINITE_VOLUME_CENTROID_COORDINATES(2:2:120);
y = x;
figure
surf(x, y, HCC_2D)
figure
surf(x, y, K1C_2D*1e-6)
figure
surf(x, y, vpv_2D)
%% close all files
fclose('all');
%% Exit message
disp('!!!!!!!!!!!!<<hardness calc End>>!!!!!!!!!!!!')

clearvars -except CHEMICAL_POTENTIALS DIMENSIONALITY DOMAIN_SIZE ...
FINITE_VOLUME_CENTROID_COORDINATES ...
GRADIENT_ENERGY_CONTRIBUTION MOLE_FRACTIONS NUMBER_OF_ELEMENTS ...
NUMBER_OF_GRID_POINTS NUMBER_OF_PHASES ...
PERMEABILITIES PHASE_FIELD PHASE_FRACTIONS TIME nel nph nts ...
ELEMENT_NEMAES PHASE_NAMES HCC work_path

disp('!!!!!!!!!!!!<<hardness calc end>>!!!!!!!!!!!!')
