
close all


%% 1D
load C:\Users\armin\Documents\Pycharmprojects\yapfi_post_processor\8-LM-2_B\vpv_liq_last_1D.TXT
load C:\Users\armin\Documents\Pycharmprojects\yapfi_post_processor\8-LM-2_B\coords_1D.TXT
LM82B_vpv = vpv_liq_last_1D;
LM82B_coords = coords_1D;
load C:\Users\armin\Documents\Pycharmprojects\yapfi_post_processor\8-LM-3_B\vpv_liq_last_1D.TXT
load C:\Users\armin\Documents\Pycharmprojects\yapfi_post_processor\8-LM-3_B\coords_1D.TXT
LM83B_vpv = vpv_liq_last_1D;
LM83B_coords = coords_1D;
load C:\Users\armin\Documents\Pycharmprojects\yapfi_post_processor\8-LM-4_B\vpv_liq_last_1D.TXT
load C:\Users\armin\Documents\Pycharmprojects\yapfi_post_processor\8-LM-4_B\coords_1D.TXT
LM84B_vpv = vpv_liq_last_1D;
LM84B_coords = coords_1D;


figure
hold on
plot(LM82B_coords, LM82B_vpv , '* red')
plot(LM83B_coords, LM83B_vpv, '+ blue')
plot(LM84B_coords, LM84B_vpv, 'o green')

%% 2D