close all;
disp('!!!!!!!!!!!!<<vpv calc Start>>!!!!!!!!!!!!')
disp(work_path)
%% General conditions Thermodynamics
N = 1; % number of moles of the system
P = 101325; % Ambient pressure
T = 1723.15; % Freezed in T
%% Dummy variables
ngp = NUMBER_OF_GRID_POINTS;
nel = NUMBER_OF_ELEMENTS;
ndim = DIMENSIONALITY;
dmsize = DOMAIN_SIZE;
mf = MOLE_FRACTIONS;
%% Open output files
if ndim > 1
    fileID_1 = fopen([work_path '\vpv_liq_first.TXT'],'w');
    fileID_2 = fopen([work_path '\vpv_liq_last.TXT'],'w');
    fileID_3 = fopen([work_path '\vpv_liq_diag_first.TXT'],'w');
    fileID_4 = fopen([work_path '\vpv_liq_line_first.TXT'],'w');
    fileID_7 = fopen([work_path '\vpv_liq_diag_last.TXT'],'w');
    fileID_8 = fopen([work_path '\vpv_liq_line_last.TXT'],'w');
    fileID_5 = fopen([work_path '\dist_diag_first.TXT'],'w');
    fileID_6 = fopen([work_path '\dist_line_first.TXT'],'w');
    fileID_9 = fopen([work_path '\dist_diag_last.TXT'],'w');
    fileID_10 = fopen([work_path '\dist_line_last.TXT'],'w');
else
    fileID_1 = fopen([work_path '\vpv_liq_first_1D.TXT'],'w');
    fileID_2 = fopen([work_path '\vpv_liq_last_1D.TXT'],'w');
    current_path= pwd;
    cd(work_path)
        copyfile FINITE_VOLUME_CENTROID_COORDINATES.TXT coords_1D.TXT
    cd(current_path)
end
%% creat a mole fraction condition string
sel =[];
for i = 1 : nel
    X_of_el(i) = cellstr([ 'X(' ELEMENT_NEMAES{i} ')']);
    if i==1
        sel = strcat(sel,ELEMENT_NEMAES{i});
    else
        sel = strcat(sel, {' '}, ELEMENT_NEMAES{i});
    end
end
%% tc calc initialization
tc_init_root;
%% Read a POLY3 file of the system
tc_read_poly3_file('wccotin_mcshp-fcca1-liquid-graphite-tcfe8.POLY3')
[nph_nw, ph_name_nw] = tc_list_phase;
[ncomp_nw, comp_name_nw] = tc_list_component;
tctoolbox('tc_set_condition', 'n', N);
tctoolbox('tc_set_condition', 'p', P);
tctoolbox('tc_set_condition', 't', T);
%% select the time step
for tstp = [1,nts]
    %% Extract mole fractions of the timestep without baoundaries
    if ndim == 1
        mf_tstp = mf((tstp-1)*nel*ngp+1: tstp*nel*ngp); 
    end
    if ndim == 2
        mf_tstp = mf((tstp-1)*nel*ngp(1)*ngp(2)+1: tstp*nel*ngp(1)*ngp(2));
    end
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
        %% read VPV of liquid 
        vpv_liq(i) = tc_get_value(char(['vpv(' char(ph_name_nw{1}) ')'])); % only Binder = liq#
        if tstp == 1
            fprintf(fileID_1, '%f \n', vpv_liq(i)); %(i) in MPa
        elseif tstp == nts
            fprintf(fileID_2, '%f \n', vpv_liq(i)); %(i) in MPa
        end
    end
%% Plotting
    %% 1D
    if ndim == 1
        x = FINITE_VOLUME_CENTROID_COORDINATES;
        figure
        plot(x, vpv_liq)
        title(['timestep' num2str(tstp)])
    end  
    %% 2D
    if ndim == 2
        %% Reshape linear arrays and surface => surfe 
        vpv_2D = reshape( vpv_liq, [ngp(1), ngp(2)] );
        x = FINITE_VOLUME_CENTROID_COORDINATES(2:2:2*ngp(1));
        y = x;
        surf(x, y, vpv_2D)
        for i = 1 : ngp(1)
            vpv_diag(i) = vpv_2D(i,i)
            xy(i) = sqrt(2*x(i)^2)
            vpv_line(i) = vpv_2D(end,i)
            %% write line and diag files
            if tstp == 1
                fprintf(fileID_3, '%f \n', vpv_diag(i)); %(i) in MPa
                fprintf(fileID_4, '%f \n', vpv_line(i)); %(i) in MPa
                fprintf(fileID_5, '%f \n', xy(i)); %(i) in MPa
                fprintf(fileID_6, '%f \n', x(i)); %(i) in MPa
            end
            if tstp == nts
                fprintf(fileID_7, '%f \n', vpv_diag(i)); %(i) in MPa
                fprintf(fileID_8, '%f \n', vpv_line(i)); %(i) in MPa
                fprintf(fileID_9, '%f \n', xy(i)); %(i) in MPa
                fprintf(fileID_10, '%f \n', x(i)); %(i) in MPa
            end
        end
        %% plot line and diag lines
        figure
        hold on
        plot(xy, vpv_diag)
        plot(x, vpv_line)
        title(['timestep' num2str(tstp)])
    end
end
%% close all files
fclose('all');
%% Exit message
disp('!!!!!!!!!!!!<<vpv_calc End>>!!!!!!!!!!!!')

clearvars -except CHEMICAL_POTENTIALS DIMENSIONALITY DOMAIN_SIZE ...
FINITE_VOLUME_CENTROID_COORDINATES ...
GRADIENT_ENERGY_CONTRIBUTION MOLE_FRACTIONS NUMBER_OF_ELEMENTS ...
NUMBER_OF_GRID_POINTS NUMBER_OF_PHASES ...
PERMEABILITIES PHASE_FIELD PHASE_FRACTIONS TIME nel nph nts ...
ELEMENT_NEMAES PHASE_NAMES HCC work_path ...
vpv_liq vpv_line xy x y ngp  vpv_2D
