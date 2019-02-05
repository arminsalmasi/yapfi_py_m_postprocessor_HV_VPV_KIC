disp('!!!!!!!!!!!!<<read profile Start>>!!!!!!!!!!!!')
disp(work_path)

%% load Yapfi text files
load ([work_path '\CHEMICAL_POTENTIALS.TXT' ],'-ascii');
load ([work_path '\DIMENSIONALITY.TXT' ],'-ascii');
load ([work_path '\DOMAIN_SIZE.TXT' ],'-ascii');
load ([work_path '\FINITE_VOLUME_CENTROID_COORDINATES.TXT' ],'-ascii');
load ([work_path '\GRADIENT_ENERGY_CONTRIBUTION.TXT' ],'-ascii');
load ([work_path '\MOLE_FRACTIONS.TXT' ],'-ascii');
load ([work_path '\NUMBER_OF_ELEMENTS.TXT' ],'-ascii');
load ([work_path '\NUMBER_OF_GRID_POINTS.TXT' ],'-ascii');
load ([work_path '\NUMBER_OF_PHASES.TXT' ],'-ascii');
load ([work_path '\PERMEABILITIES.TXT' ],'-ascii');
load ([work_path '\PHASE_FIELD.TXT' ],'-ascii');
load ([work_path '\PHASE_FRACTIONS.TXT' ],'-ascii');
load ([work_path '\TIME.TXT' ],'-ascii');
nts = size(TIME,1);
nel = NUMBER_OF_ELEMENTS;
nph = NUMBER_OF_PHASES;
%% Read element naems word buy word
[fid,msg] = fopen([work_path '\ELEMENT_NAMES.TXT' ],'r');
[val,count] = fread(fid);
res = fclose(fid);
ELEMENT_NEMAES = cell(0);
j=1;
for i = 1: nel
  k = j;
  ok = 1;
  while ok==1
    k = k+1;
    if val(k)==13
      val(k) = [];
    end
    if val(k)==10
      ok = 0;
    end
  end
  ELEMENT_NEMAES{i} = char(val(j:k-1)');
  j=k+1;
end
%% read phase namEs word by word
[fid,msg] = fopen([work_path '\PHASE_NAMES.TXT' ],'r');
[val,count] = fread(fid);
res = fclose(fid);
PHASE_NAMES = cell(0);
j = 1;
for i = 1: nph
  k = j; 
  ok = 1;
  while ok==1
    k = k+1;
    if val(k)==13
      val(k) = [];
    end
    if val(k)==10
      ok = 0;
    end
  end
  PHASE_NAMES{i} = char(val(j:k-1)');
  j = k+1;
end
%% clear unnecessary variables
clearvars -except CHEMICAL_POTENTIALS DIMENSIONALITY DOMAIN_SIZE ...
FINITE_VOLUME_CENTROID_COORDINATES ...
GRADIENT_ENERGY_CONTRIBUTION MOLE_FRACTIONS NUMBER_OF_ELEMENTS ...
NUMBER_OF_GRID_POINTS NUMBER_OF_PHASES ...
PERMEABILITIES PHASE_FIELD PHASE_FRACTIONS TIME nel nph nts ...
ELEMENT_NEMAES PHASE_NAMES work_path
%% Exit message

disp('!!!!!!!!!!!!<<read profile end>>!!!!!!!!!!!!')

%% next script
%hardness_WTICNCo
vpv_calc
%dictra_test
