disp(work_path)

fileID_1 = fopen([work_path '\timestep_molefraction.TXT'],'w');

ngp = NUMBER_OF_GRID_POINTS;
nel = NUMBER_OF_ELEMENTS;
ndim = DIMENSIONALITY;
dmsize = DOMAIN_SIZE;
mf = MOLE_FRACTIONS;

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
T = 1723.15; % Freezed in T
%%
tc_init_root;
tc_read_poly3_file('wccotin_mcshp-fcca1-tcfe8.POLY3')
%% select the time step
tstp = 1;
%% Extract mole fractions of the timestep without baoundaries
%% this will be taken care of in the python postprocessor
if ndim == 1
    mf_tstp = mf((tstp-1)*nel*ngp(1)+1: tstp*nel*ngp(1));
end
fprintf(fileID_1, '%f \n', mf_tstp(:))

dic_command('go -m')
dic_command('s-cond glob t 0 1723.15; * N')
dic_command('en-reg cc')
dic_command('en-gr cc 1e-3 r-p-b-p 70 FINITE_VOLUME_CENTROID_COORDINATES.DAT')
dic_command('en-ph ACTIVE CC MATRIX liq')





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
