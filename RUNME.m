MIMM_set_path 

%% Dictionary Generation
%This step takes approximately 1-2 days. I suggest to use of the provided
%precalculated dictionaries.

% Stochastic Dictionary
dictionary = generate_dictionary_stochastic();
save("Dictionary\MIMM_dictionary_stochastic.mat","dictionary")

% Deterministic Dictionary
dictionary = generate_dictionary_deterministic();
save("Dictionary\MIMM_dictionary_deterministic.mat","dictionary")



%% Matching Pursuit
load('FA.mat')
load('QSM.mat')
load('theta.mat')
load('iField.mat','iField','TE')


%% Load the stochastic dictionary
load('MIMM_dictionary_stochastic.mat')

% Basic MIMM
orientation_strategy = "basic"; 
MIMM_basic = MIMM(dictionary, QSM, Brain_Mask, iField, TE,orientation_strategy);

% DTI Orientation Informed MIMM
orientation_strategy = "orientation_informed"; 
MIMM_DTI = MIMM(dictionary, QSM, Brain_Mask, iField, TE,orientation_strategy,FA_DTI,theta_DTI);

% Atlas Orientation Informed MIMM

orientation_strategy = "orientation_informed"; 
MIMM_Atlas = MIMM(dictionary, QSM, Brain_Mask, iField, TE,orientation_strategy,FA_atlas,theta_atlas);

save("Example_Results\stochastic_MIMM_results.mat","MIMM_Atlas","MIMM_DTI","MIMM_basic")


%% Load the deterministic dictionary

load('MIMM_dictionary_deterministic.mat')

% Basic MIMM
orientation_strategy = "basic"; 
MIMM_basic = MIMM(dictionary, QSM, Brain_Mask, iField, TE,orientation_strategy);

% DTI Orientation Informed MIMM
orientation_strategy = "orientation_informed"; 
MIMM_DTI = MIMM(dictionary, QSM, Brain_Mask, iField, TE,orientation_strategy,FA_DTI,theta_DTI);

% Atlas Orientation Informed MIMM

orientation_strategy = "orientation_informed"; 
MIMM_Atlas = MIMM(dictionary, QSM, Brain_Mask, iField, TE,orientation_strategy,FA_atlas,theta_atlas);

save("Example_Results\deterministic_MIMM_results.mat","MIMM_Atlas","MIMM_DTI","MIMM_basic")