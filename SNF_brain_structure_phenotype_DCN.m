%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This script contains the code for the paper,
%   Hong and Sisk et al, 2021 "Decomposing complex links between the childhood environment and brain structure in school-aged youth". 
%
% - It generates necessary components to generate the figures. Before running, multiple matlab toolboxes should be installed 
%   (the links are appended below). 
% - Because the entire code is quite long, if you want to see it with some structures, please fold the codes first.
%   (https://www.mathworks.com/help/matlab/matlab_prog/improve-code-readability.html#brazeg3-1)
%
% - Last updated: Feb 18, 2021 (hong.seok.jun@gmail.com)


%% 1) Determine dataset
clear all
close all

addpath(genpath('surfstat'));               % https://mica-mni.github.io/surfstat/
addpath(genpath('SNFmatlab'));              % http://compbio.cs.toronto.edu/SNF/SNF/Software_files/SNFmatlab_v2.1.zip
addpath(genpath('ComBatHarmonization'));    % https://github.com/Jfortin1/ComBatHarmonization
addpath(genpath('export_fig'));             % https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig
addpath(genpath('ComDetTBv090'));           % https://www.mathworks.com/matlabcentral/fileexchange/45867-community-detection-toolbox
addpath(genpath('gifti-master'));           % https://github.com/gllmflndn/gifti/archive/master.zip
addpath(genpath('BCT'));                    % https://sites.google.com/site/bctnet/Home/functions/BCT.zip?attredirects=0
addpath(genpath('effect_size'));            % https://www.mathworks.com/login?uri=https%3A%2F%2Fwww.mathworks.com%2Fmatlabcentral%2Ffileexchange%2F32398-hhentschke-measures-of-effect-size-toolbox&form_type=community
addpath(genpath('Violinplot-Matlab'));      % https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot

for colormap_preparation = 1
    
    mycol.blackblue =   flipud( [zeros(1,3)*0.8; zeros(127,1) ...
        (0:126)'/127 ones(127,1)]);
    mycol.blue      =   flipud( [ones(1,3)*0.8; ...
        zeros(127,1) (0:126)'/127 ones(127,1)]);
    mycol.red       =   flipud( [ones(1,3)*0.8; ...
        ones(500,1) linspace(0,253,500)'/254 zeros(500,1);...
        ones(64,1) ones(64,1) linspace(0,200,64)'/254]);
    
    yeo_colormap = [ 0 0 0;
        120 18 134;
        70 130 180;
        0 118 14;
        196 58 250;
        220 248 164;
        230 148 34;
        205 62 78 ]/255;
    
end
   
for data_load_clean = 1
    
    % The sources of the phenotypic data are complex because of their history in building the excel file.
    % There are three difference sources.
    % 1) The original database from what I and Anthony have worked together. This is a split version for 1.1 and 2.0. These are snf_discovery.mat and snf_replication.mat.
    % 2) The next one is the one that Lucinda has created. The name of files starts with LS_PhenoData_XX. 
    %    There are three files (one for all, the second for summed trauma scores and the final one for summed trauma and family history). 
    % 3) The third one is what Anthony has created by the request of Lucinda for including race and ethnic information.
    % 
    % Because of these multi-sources of data, there needs to be a merging process with subject ID as an anchor point. 
    snf_discovery   = load('snf_discovery.mat');
    snf_replication = load('snf_replication.mat');
    LSPhenoDataSample_disc_repl = importfile_phenoData_TraumaSum_FamHist_allVars([ 'LS_PhenoData_TraumaSum_FamHist_allIDs_2.29.20.csv' ], 2, 1001); clinical_idx = 1:8; validation_idx = 9:11;
    fulldemographics1000fromanthony = importfile_fulldemo('full_demographics_1000_from_anthony.csv', 2, 1001);
    LShouseincome_disc_repl = importfile_houseincome('ABCD_family_income.csv');
    cbclscases = importfile_cbcls_cases("/Users/hongsukjun/Documents/03_data/01_documents/02_research/00_Project/17_ABCD_Environmental/SNF/codes/cbcls_cases_092120.csv", [2, Inf]);
    ST = snf_discovery.ST; 
    mask = snf_discovery.mask;
    
    yeo_krienen_map_left  = SurfStatReadData('Yeo_parcellation_conte69mapping_left.txt');
    yeo_krienen_map_right = SurfStatReadData('Yeo_parcellation_conte69mapping_right.txt');
    yeo_krienen_map = [yeo_krienen_map_left yeo_krienen_map_right]; yeo_krienen_map(~mask) = 0;
    
    % 1-1) snf vs. LSPheno matching
    CASE = [ snf_discovery.CASE; snf_replication.CASE; ];
    CASE_match = [];
    for i = 1 : size(CASE)
        
       CASE_match = [ CASE_match; find(strcmp(LSPhenoDataSample_disc_repl.subjectkey, CASE{i})) ];
        
    end
    
    % 1-2) sanity check: should be zero.
    crpbi_y_ss_parent1 = [ snf_discovery.crpbi_y_ss_parent1; snf_replication.crpbi_y_ss_parent1; ];
    sum(crpbi_y_ss_parent1-LSPhenoDataSample_disc_repl.crpbi_y_ss_parent(CASE_match))
    clear('crpbi_y_ss_parent1');
    
    % 1-3) phenotypic score and matching and build-up 
    %      subtyping variables: 
    %         crpbi_y_ss_parent         - Children's Report of Parental Behavioral Inventory (Acceptance Subscale Mean of Report by Parent): 1) Makes me feel better, 2) Smiles at me very often, 3) Is able to make me feel better when I am upset...
    %         crpbi_y_ss_caregiver      - Children's Report of Parental Behavioral Inventory (Acceptance Subscale Mean of Report by Secondary Caregiver): 1) Makes me feel better, 2) Smiles at me very often, 3) Is able to make me feel better when I am upset...
    %         physical_abuse_sum        - Physical abuse total
    %         sexual_abuse_sum          - Sexual abuse total
    %         exposure_violence_sum     - Violence exposure total
    %         trauma_other_sum          - Other trauma total
    %         mentalhealth_sum          - Mental health total
    %         substance_sum             - Substance use total
    %         criminal_sum              - Criminal total
    %         pmq_y_ss_mean             - Parental Monitoring: 1) How often do your parents/guardians know where you are? 2) How often do your parents know who you are with when you are not at school and away from home? and so on.
    %         srpf_y_ss_ses             - School Environment Subscale (e.g. In my school, students have lots of chances to help decide things like class activities and rules, 	I get along with my teachers)
    %         nsc_p_ss_mean_3_items     - Neighborhood Safety Protocol
    %         fes_p_ss_fc               - Conflict subscale from the Family Environment Scale Sum of Parent Report 
    %
    %      validation variables:
    %         cbcl_scr_syn_internal_t   - Child Behavior Checklist Scores: Internal CBCL Syndrome Scale (t-score)
    %         cbcl_scr_syn_external_t   - Child Behavior Checklist Scores: External CBCL Syndrome Scale (t-score)
    %         cbcl_scr_syn_totprob_t    - Child Behavior Checklist Scores: TotProb CBCL Syndrome Scale (t-score)
    fieldnames_set = LSPhenoDataSample_disc_repl.Properties.VariableNames;
    phenotypic_data = [];
    for i = 1 : length(fieldnames_set)
        
        eval([ fieldnames_set{i} ' = LSPhenoDataSample_disc_repl.' fieldnames_set{i} '(CASE_match);' ]);
        eval([ fieldnames_set{i} '_discovery   = ' fieldnames_set{i} '(randidx(1:500));' ]);
        eval([ fieldnames_set{i} '_replication = ' fieldnames_set{i} '(randidx(501:1000));' ]);
        
        if(i>1)
            eval([ 'phenotypic_data = [ phenotypic_data ' fieldnames_set{i} '];' ]);
        end
        
    end
    
    phenotypic_data               = phenotypic_data(randidx, :);
    phenotypic_data_discovery     = phenotypic_data(1:500, :);
    phenotypic_data_replication   = phenotypic_data(501:1000, :);
    cbcl_subsale_data             = cbclscases(randidx, :);
    cbcl_subsale_data_discovery   = cbcl_subsale_data(1:500, :);
    cbcl_subsale_data_replication = cbcl_subsale_data(501:1000, :);
    
    % 1-4) snf vs. fulldemoAnthony matching
    CASE_match = [];
    SUBJECTKEY = erase(fulldemographics1000fromanthony.SUBJECTKEY, '_');
    for i = 1 : size(CASE)
        
       CASE_match = [ CASE_match; find(strcmp(SUBJECTKEY, CASE{i})) ];
        
    end
    clear('SUBJECTKEY');
    
    % 1-5) sanity check: should be zero.
    AGE = [ snf_discovery.AGE; snf_replication.AGE; ];
    sum(AGE-fulldemographics1000fromanthony.INTERVIEW_AGE(CASE_match))
    clear('AGE');
    
    % 1-6) demographic score and matching and build-up 
    %         SEX               
    %         INTERVIEW_AGE     
    %         SITE_ID_L         
    %         Release           
    %         DEMO_RACE_A_P___10 - White Blanca (child)             1
    %         DEMO_RACE_A_P___11 - Black/African (child)            2
    %         DEMO_RACE_A_P___12 - American Indian (child)          3
    %         DEMO_RACE_A_P___13 - Alaska Native (child)            4
    %         DEMO_RACE_A_P___14 - Native Hawaiian (child)          5
    %         DEMO_RACE_A_P___15 - Guamanian (child)                6
    %         DEMO_RACE_A_P___16 - Samoan (child)                   7
    %         DEMO_RACE_A_P___17 - Other Pacific Islander (child)   8
    %         DEMO_RACE_A_P___18 - Asian Indian (child)             9
    %         DEMO_RACE_A_P___19 - Chinese (child)                  10
    %         DEMO_RACE_A_P___20 - Filipino (child)                 11
    %         DEMO_RACE_A_P___21 - Japanese (child)                 12
    %         DEMO_RACE_A_P___22 - Korean (child)                   13
    %         DEMO_RACE_A_P___23 - Vietnamese (child)               14
    %         DEMO_RACE_A_P___24 - Other Asian (child)              15
    %         DEMO_RACE_A_P___25 - Other Race (child)               16
    %         DEMO_RACE_A_P___77 - Refuse To Answer                 17
    %         DEMO_RACE_A_P___99 - Don't Know                       18
    %         DEMO_ETHN_V2       - Do you consider the child Hispanic/Latino/Latina?   
    %         DEMO_PRNT_ED_V2    - What is the highest grade or level of school you have completed or the highest degree you have received? 
    %         DEMO_PRNT_PRTNR_V2 - Do you have a partner?
    %         DEMO_PRTNR_ED_V2   - What is the highest grade or level of school your partner completed or highest degree they received?
    %         DEMO_PRNT_INCOME_V2  - How much did you earn, before taxes and other deductions, during the past 12 months?
    %         DEMO_PRTNR_INCOME_V2 - How much does your partner earn, before taxes and other deductions, during the past 12 months?
    fieldnames_set = fulldemographics1000fromanthony.Properties.VariableNames;
    for i = 1 : length(fieldnames_set)
        
        eval([ fieldnames_set{i} ' = fulldemographics1000fromanthony.' fieldnames_set{i} '(CASE_match);' ]);
        eval([ fieldnames_set{i} '_discovery   = ' fieldnames_set{i} '(randidx(1:500));' ]);
        eval([ fieldnames_set{i} '_replication = ' fieldnames_set{i} '(randidx(501:1000));' ]);
        
    end
    
    demorace_data = [ DEMO_RACE_A_P___10 DEMO_RACE_A_P___11 DEMO_RACE_A_P___12 DEMO_RACE_A_P___13 DEMO_RACE_A_P___14 DEMO_RACE_A_P___15 ...
                      DEMO_RACE_A_P___16 DEMO_RACE_A_P___17 DEMO_RACE_A_P___18 DEMO_RACE_A_P___19 DEMO_RACE_A_P___20 ...
                      DEMO_RACE_A_P___21 DEMO_RACE_A_P___22 DEMO_RACE_A_P___23 DEMO_RACE_A_P___24 DEMO_RACE_A_P___25 ...
                      DEMO_RACE_A_P___77 DEMO_RACE_A_P___99 ];
    demorace_name = { 'WhiteBlanca', 'BlackAfrican', 'AmericanIndian', 'AlaskaNative', 'NativeHawaiian', 'Guamanian', 'Samoan', ...
                      'OtherPacificIslander', 'AsianIndian', 'Chinese', 'Filipino', 'Japanese', 'Korean', 'Vietnamese', 'OtherAsian', ...
                      'OtherRace', 'RefuseToAnswer', 'DonotKnow' };
    demorace_data = demorace_data(randidx, :);
    
    houseincome_data = LShouseincome_disc_repl.demo_comb_income_v2;
    
    count = 1;
    demorace_lookup_name = cell(1, 1);
    demorace_lookup      = [];
    demorace_categorical = zeros(size(demorace_data, 1), 1);
    for i = 1 : size(demorace_data, 1)
        
        if(demorace_categorical(i) == 0)
            
            demorace_lookup_name{count, 1} =  strjoin(demorace_name(logical(demorace_data(i, :))), ', ');
            demorace_lookup = [ demorace_lookup; count demorace_data(i, :) ];
            demorace_categorical(find(sum((demorace_data - repmat(demorace_data(i, :), size(demorace_data, 1), 1)).^2, 2) == 0)) = count;
            count = count + 1;
            
        end
        
    end
    
    demographic_data = [ SEX=='F' INTERVIEW_AGE SITE_ID_L demorace_categorical DEMO_ETHN_V2 ...
                         DEMO_PRNT_ED_V2 DEMO_PRNT_PRTNR_V2 DEMO_PRNT_INCOME_V2 DEMO_PRNT_INCOME_V2 ];
    
    demographic_data = demographic_data(randidx, :);                 
    demographic_data_discovery = demographic_data(1:500, :);
    demographic_data_replication = demographic_data(501:1000, :);
    
    DEMO_PRNT_INCOME_V2(isnan(DEMO_PRNT_INCOME_V2)) = 999;
    DEMO_PRTNR_INCOME_V2(isnan(DEMO_PRTNR_INCOME_V2)) = 999;
    for i = 1 : size(DEMO_PRNT_INCOME_V2, 1)
        if(DEMO_PRNT_INCOME_V2(i) > 21 & DEMO_PRTNR_INCOME_V2(i) > 21)
            
            DEMO_HOUSE_INCOME(i) = NaN;
            
        elseif(DEMO_PRNT_INCOME_V2(i) > 21 & DEMO_PRTNR_INCOME_V2(i) <= 21)
            
            DEMO_HOUSE_INCOME(i) = DEMO_PRTNR_INCOME_V2(i);
            
        elseif(DEMO_PRNT_INCOME_V2(i) <= 21 & DEMO_PRTNR_INCOME_V2(i) > 21)
            
            DEMO_HOUSE_INCOME(i) = DEMO_PRNT_INCOME_V2(i);
            
        else
            
            DEMO_HOUSE_INCOME(i) = max(DEMO_PRNT_INCOME_V2(i), DEMO_PRTNR_INCOME_V2(i));
            
        end
        
    end
    
    probability_dist = [];
    for i = 1 : size(demographic_data_discovery, 2)
        
        freq = tabulate(demographic_data_discovery(:, i));
        prob = 1./(freq(:, 3)/100); prob = prob / sum(prob(~isinf(prob)));
        prob_ind = zeros(size(demographic_data_discovery, 1), 1);
        for j = 1 : size(freq, 1)
           
            prob_ind(demographic_data_discovery(:, i) == freq(j, 1)) = prob(j);
            
        end
        probability_dist = [ probability_dist prob_ind ];
        
    end
    
    fn_prob_discovery = probability_dist(:, 1);
    for i = 2 : size(probability_dist, 2)
       
        fn_prob_discovery = fn_prob_discovery.*probability_dist(:, i);
        
    end
    fn_prob_discovery = fn_prob_discovery / sum(fn_prob_discovery);
    
    probability_dist = [];
    for i = 1 : size(demographic_data_replication, 2)
        
        freq = tabulate(demographic_data_replication(:, i));
        prob = 1./(freq(:, 3)/100); prob = prob / sum(prob(~isinf(prob)));
        prob_ind = zeros(size(demographic_data_replication, 1), 1);
        for j = 1 : size(freq, 1)
           
            prob_ind(demographic_data_replication(:, i) == freq(j, 1)) = prob(j);
            
        end
        probability_dist = [ probability_dist prob_ind ];
        
    end
    
    fn_prob_replication = probability_dist(:, 1);
    for i = 2 : size(probability_dist, 2)
       
        fn_prob_replication = fn_prob_replication.*probability_dist(:, i);
        
    end
    fn_prob_replication = fn_prob_replication / sum(fn_prob_replication);
    
    cbclscore_name = { 'CBCL_SCR_SYN_ANXDEP_T', 'CBCL_SCR_SYN_WITHDEP_T', 'CBCL_SCR_SYN_SOMATIC_T', 'CBCL_SCR_SYN_SOCIAL_T', 'CBCL_SCR_SYN_THOUGHT_T', ...
                       'CBCL_SCR_SYN_ATTENTION_T', 'CBCL_SCR_SYN_RULEBREAK_T', 'CBCL_SCR_SYN_AGGRESSIVE_T', 'CBCL_SCR_DSM5_DEPRESS_T', 'CBCL_SCR_DSM5_ANXDISORD_T', ...
                       'CBCL_SCR_DSM5_SOMATICPR_T', 'CBCL_SCR_DSM5_ADHD_T', 'CBCL_SCR_DSM5_OPPOSIT_T', 'CBCL_SCR_DSM5_CONDUCT_T', 'CBCL_SCR_07_SCT_T', ' CBCL_SCR_07_OCD_T', 'CBCL_SCR_07_STRESS_T'};
                   
    cbcl_subscale = [];
    for i = 1 : size(cbclscore_name, 2)
        
        eval([ 'cbcl_subscale = [ cbcl_subscale cbcl_subsale_data.' cbclscore_name{i} ' ];' ]);
        
    end
    cbcl_subscale_disc = cbcl_subscale(1:500, :);
    cbcl_subscale_repl = cbcl_subscale(501:1000, :);
    
end

%% 2) Regress out aging, sex and site effects using linear model and ComBat approaches
for confounding_factor_statistical_correction = 1
    
    thickness_set_org = [ snf_discovery.thickness_set_org; snf_replication.thickness_set_org ]; thickness_set_org = thickness_set_org(randidx, :);
    myelin_set_org = [ snf_discovery.myelin_set_org; snf_replication.myelin_set_org ]; myelin_set_org = myelin_set_org(randidx, :);
    CASE_org = [ snf_discovery.CASE; snf_replication.CASE ]; CASE_org = CASE_org(randidx);
    AGE_org = [ snf_discovery.AGE; snf_replication.AGE ]; AGE_org = AGE_org(randidx);
    GENDER_org = [ snf_discovery.GENDER; snf_replication.GENDER ]; GENDER_org = GENDER_org(randidx);
    SITE_org = [ snf_discovery.SITE; snf_replication.SITE ]; SITE_org = SITE_org(randidx);
    demorace_data_org = demorace_data; demorace_data_org = demorace_data_org(randidx, :);
    ethn_data_org = DEMO_ETHN_V2; ethn_data_org = ethn_data_org(randidx);
    parental_edu_org = DEMO_PRNT_ED_V2; parental_edu_org = parental_edu_org(randidx);
    parental_income_org = DEMO_HOUSE_INCOME; parental_income_org = parental_income_org(randidx);
    houseincome_org = houseincome_data; houseincome_org = houseincome_org(randidx);
    
    % 2-1) discovery process
    for discovery_set = 1
        
        % 2-1-1) variable setup
        fieldnames_set             = LSPhenoDataSample_disc_repl.Properties.VariableNames;

        thickness_set_org_disc     = thickness_set_org(1:500, :);
        myelin_set_org_disc        = myelin_set_org(1:500, :);
        clinicalvar_set_org_disc   = phenotypic_data_discovery(:, clinical_idx);
        clinicalvar_names_disc     = fieldnames_set(clinical_idx+1);
        clinicalvar2_set_org_disc  = cbcl_subscale_disc;
        clinicalvar2_names_disc    = cbclscore_name;
        validationvar_set_org_disc = phenotypic_data_discovery(:, validation_idx);
        validationvar_names_disc   = fieldnames_set(validation_idx+1);
        demorace_data_org_disc     = demorace_data_org(1:500, :);
        ethn_data_org_disc         = ethn_data_org(1:500, :);
        parental_edu_org_disc      = parental_edu_org(1:500, :);
        parental_income_org_disc   = parental_income_org(1:500);
        houseincome_org_disc       = houseincome_org(1:500); 
        
        CASE_FN_disc       = CASE_org(1:500);
        AGE_FN_disc        = AGE_org(1:500);
        GENDER_FN_disc     = GENDER_org(1:500);
        SITE_FN_disc       = cellstr(num2str(SITE_org(1:500)));

        thickness_set_disc       = thickness_set_org_disc;
        myelin_set_disc          = myelin_set_org_disc;
        clinicalvar_set_disc     = clinicalvar_set_org_disc;
        validationvar_set_disc   = validationvar_set_org_disc;
        clinicalvar2_set_disc    = clinicalvar2_set_org_disc;
        
        demorace_data_set_disc   = demorace_data_org_disc;
        ethn_data_set_disc       = ethn_data_org_disc;
        parental_edu_set_disc    = parental_edu_org_disc;
        parental_income_set_disc = parental_income_org_disc; 
        houseincome_set_disc     = houseincome_org_disc;

        pheno_scores_disc = [ clinicalvar_set_disc validationvar_set_disc ];
        temp = corr(pheno_scores_disc);
        
        figure; heatmap(temp, fieldnames_set([clinical_idx validation_idx]+1), fieldnames_set([clinical_idx validation_idx]+1)); xtickangle(45); colormap('jet'); caxis([-0.5 0.5]);
        
        % 2-1-2) outlier check
        threshold_for_outlier = 0.2;
        
        % - imaging data
        zthickness_set_org_disc = (thickness_set_org_disc - repmat(mean(thickness_set_org_disc, 1), size(thickness_set_org_disc, 1), 1))./repmat(std(thickness_set_org_disc, 0, 1), size(thickness_set_org_disc, 1), 1);
        thresholded_thickness_disc = abs(zthickness_set_org_disc)>3.09;
        thresholded_thickness_disc = sum(thresholded_thickness_disc, 2)/size(thresholded_thickness_disc, 2);
        [a b] = sort(thresholded_thickness_disc); a'
        outlier_idx_disc = unique(b(a>=threshold_for_outlier));
        
        zmyelin_set_org_disc = (myelin_set_org_disc - repmat(mean(myelin_set_org_disc, 1), size(myelin_set_org_disc, 1), 1))./repmat(std(myelin_set_org_disc, 0, 1), size(myelin_set_org_disc, 1), 1);
        thresholded_myelin_disc = abs(zmyelin_set_org_disc)>3.09;
        thresholded_myelin_disc = sum(thresholded_myelin_disc, 2)/size(thresholded_myelin_disc, 2);
        [a b] = sort(thresholded_myelin_disc); a'
        outlier_idx_disc = unique([ outlier_idx_disc; b(a>=threshold_for_outlier) ]);
        
        thickness_set_disc(outlier_idx_disc, :)     = [];
        myelin_set_disc(outlier_idx_disc, :)        = [];
        clinicalvar_set_disc(outlier_idx_disc, :)   = [];
        validationvar_set_disc(outlier_idx_disc, :) = [];
        clinicalvar2_set_disc(outlier_idx_disc, :)   = [];
        CASE_FN_disc(outlier_idx_disc)              = [];
        AGE_FN_disc(outlier_idx_disc)               = [];
        GENDER_FN_disc(outlier_idx_disc)            = [];
        SITE_FN_disc(outlier_idx_disc)              = [];
        fn_prob_discovery(outlier_idx_disc)         = [];
        demographic_data_discovery(outlier_idx_disc, :) = [];
        demorace_data_set_disc(outlier_idx_disc, :) = [];
        ethn_data_set_disc(outlier_idx_disc)        = [];
        parental_edu_set_disc(outlier_idx_disc)     = [];
        parental_income_set_disc(outlier_idx_disc)     = [];
        Release_discovery(outlier_idx_disc)         = [];
        houseincome_set_disc(outlier_idx_disc)         = [];
        
        % 2-1-3) statistically remove the age, sex and site effects for brain imaging data and the site effect for phenotype scores
        for age_sex_correction = 1
            
            for thickness = 1
                
                AGE_term = term(AGE_FN_disc);
                SEX_term = term(GENDER_FN_disc);
                M = 1 + AGE_term + SEX_term;
                slm  = SurfStatLinMod(thickness_set_disc, M, ST);
                %                 slm  = SurfStatT(slm, AGE_FN_disc);
                %                 slm  = SurfStatT(slm, SEX_term.M-SEX_term.F);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(mean(thickness_set_disc, 1), ST);  SurfStatColLim([0 4.5]);
                
                thickness_set_disc_age_sex_corr = thickness_set_disc - slm.X(:, 2:4)*slm.coef(2:4, :);
                slm  = SurfStatLinMod(thickness_set_disc_age_sex_corr, M, ST);
                %                 slm  = SurfStatT(slm, AGE_FN_disc);
                %                 slm  = SurfStatT(slm, SEX_term.M-SEX_term.F);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(mean(thickness_set_disc_age_sex_corr, 1), ST);  SurfStatColLim([0 4.5]);
                
            end
            for myelin = 1
                
                AGE_term = term(AGE_FN_disc);
                SEX_term = term(GENDER_FN_disc);
                M = 1 + AGE_term + SEX_term;
                slm  = SurfStatLinMod(myelin_set_disc, M, ST);
                %                 slm  = SurfStatT(slm, AGE_FN_disc);
                %                 slm  = SurfStatT(slm, SEX_term.M-SEX_term.F);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(mean(myelin_set_disc, 1), ST);  SurfStatColLim([0.3 0.6]);
                
                myelin_set_disc_age_sex_corr = myelin_set_disc - slm.X(:, 2:4)*slm.coef(2:4, :);
                slm  = SurfStatLinMod(myelin_set_disc_age_sex_corr, M, ST);
                %                 slm  = SurfStatT(slm, AGE_FN_disc);
                %                 slm  = SurfStatT(slm, SEX_term.M-SEX_term.F);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(mean(myelin_set_disc_age_sex_corr, 1), ST); SurfStatColLim([0.25 0.35]);
                
            end
            
        end
        for site_correction = 1
            
            for thickness = 1
                
                thickness_set_std = std(thickness_set_disc_age_sex_corr, 0, 1);
                thickness_set_temp = thickness_set_disc_age_sex_corr;
                
                thickness_set_temp(:, find(thickness_set_std == 0)) = [];
                
                mod = [];
                bayes_thickness_set_temp = combat(thickness_set_temp', str2num(cell2mat(SITE_FN_disc)), mod, 1)';
                bayes_thickness_set_disc = zeros(length(CASE_FN_disc), 20484); bayes_thickness_set_disc(:,  find(thickness_set_std ~= 0)) = bayes_thickness_set_temp;
                
                SITE_term = term(SITE_FN_disc);
                M1 = 1 + SITE_term;
                M2 = 1;
                slm1  = SurfStatLinMod(thickness_set_disc_age_sex_corr, M1, ST);
                slm2  = SurfStatLinMod(thickness_set_disc_age_sex_corr, M2, ST);
                %                 slm  = SurfStatF(slm1, slm2);
                %                 [ pval, peak, clus, clusid ] = SurfStatP(slm, mask, 0.025);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(pval.C, ST); SurfStatColLim([0 0.05]);
                
                SITE_term = term(SITE_FN_disc);
                M1 = 1 + SITE_term;
                M2 = 1;
                slm1  = SurfStatLinMod(bayes_thickness_set_disc, M1, ST);
                slm2  = SurfStatLinMod(bayes_thickness_set_disc, M2, ST);
                %                 slm  = SurfStatF(slm1, slm2);
                %                 [ pval, peak, clus, clusid ] = SurfStatP(slm, mask, 0.025);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(pval.C, ST); SurfStatColLim([0 0.05]);
                
            end
            for myelin = 1
                
                myelin_set_std = std(myelin_set_disc_age_sex_corr, 0, 1);
                myelin_set_temp = myelin_set_disc_age_sex_corr;
                
                myelin_set_temp(:, find(myelin_set_std == 0)) = [];
                
                mod = [];
                bayes_myelin_set_temp = combat(myelin_set_temp', str2num(cell2mat(SITE_FN_disc)), mod, 1)';
                bayes_myelin_set_disc = zeros(length(CASE_FN_disc), 20484); bayes_myelin_set_disc(:,  find(myelin_set_std ~= 0)) = bayes_myelin_set_temp;
                
                SITE_term = term(SITE_FN_disc);
                M1 = 1 + SITE_term;
                M2 = 1;
                slm1  = SurfStatLinMod(myelin_set_disc_age_sex_corr, M1, ST);
                slm2  = SurfStatLinMod(myelin_set_disc_age_sex_corr, M2, ST);
                %                 slm  = SurfStatF(slm1, slm2);
                %                 [ pval, peak, clus, clusid ] = SurfStatP(slm, mask, 0.025);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(pval.C, ST); SurfStatColLim([0 0.05]);
                
                SITE_term = term(SITE_FN_disc);
                M1 = 1 + SITE_term;
                M2 = 1;
                slm1  = SurfStatLinMod(bayes_myelin_set_disc, M1, ST);
                slm2  = SurfStatLinMod(bayes_myelin_set_disc, M2, ST);
                %                 slm  = SurfStatF(slm1, slm2);
                %                 [ pval, peak, clus, clusid ] = SurfStatP(slm, mask, 0.025);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(pval.C, ST); SurfStatColLim([0 0.05]);
                
            end
            for clinicalvar = 1
                
                clinicalvar_set_std   = std(clinicalvar_set_disc, 0, 1);
                clinicalvar_set_temp  = clinicalvar_set_disc;
                clinicalvar_set_temp(:, find(clinicalvar_set_std == 0)) = [];
                clinicalvar2_set_std  = std(clinicalvar2_set_disc, 0, 1);
                clinicalvar2_set_temp = clinicalvar2_set_disc;
                clinicalvar2_set_temp(:, find(clinicalvar2_set_std == 0)) = [];
                
                mod = [];
                bayes_clinicalvar_set_temp = combat(clinicalvar_set_temp', str2num(cell2mat(SITE_FN_disc)), mod, 1)';
                bayes_clinicalvar_set_disc = zeros(length(CASE_FN_disc), size(clinicalvar_set_disc, 2)); bayes_clinicalvar_set_disc(:,  find(clinicalvar_set_std ~= 0)) = bayes_clinicalvar_set_temp;
                bayes_clinicalvar2_set_temp = combat(clinicalvar2_set_temp', str2num(cell2mat(SITE_FN_disc)), mod, 1)';
                bayes_clinicalvar2_set_disc = zeros(length(CASE_FN_disc), size(clinicalvar2_set_disc, 2)); bayes_clinicalvar2_set_disc(:,  find(clinicalvar2_set_std ~= 0)) = bayes_clinicalvar2_set_temp;
                %                 SITE_term = term(SITE_FN_disc);
                %                 M1 = 1 + SITE_term;
                %                 M2 = 1;
                %                 slm1  = SurfStatLinMod(clinicalvar_set_disc, M1);
                %                 slm2  = SurfStatLinMod(clinicalvar_set_disc, M2);
                %                 slm  = SurfStatF(slm1, slm2); slm.t

                %                 SITE_term = term(SITE_FN_disc);
                %                 M1 = 1 + SITE_term;
                %                 M2 = 1;
                %                 slm1  = SurfStatLinMod(bayes_clinicalvar_set_disc, M1);
                %                 slm2  = SurfStatLinMod(bayes_clinicalvar_set_disc, M2);
                %                 slm  = SurfStatF(slm1, slm2); slm.t
                
            end
            for validationlvar = 1
                
                validationvar_set_std = std(validationvar_set_disc, 0, 1);
                validationvar_set_temp = validationvar_set_disc;
                validationvar_set_temp(:, find(validationvar_set_std == 0)) = [];
                
                mod = [];
                bayes_validationvar_set_temp = combat(validationvar_set_temp', str2num(cell2mat(SITE_FN_disc)), mod, 1)';
                bayes_validationvar_set_disc = zeros(length(str2num(cell2mat(SITE_FN_disc))), size(validationvar_set_disc, 2)); bayes_validationvar_set_disc(:,  find(validationvar_set_std ~= 0)) = bayes_validationvar_set_temp;
                %                 SITE_term = term(SITE_FN_disc);
                %                 M1 = 1 + SITE_term;
                %                 M2 = 1;
                %                 slm1  = SurfStatLinMod(validationvar_set_disc, M1);
                %                 slm2  = SurfStatLinMod(validationvar_set_disc, M2);
                %                 slm  = SurfStatF(slm1, slm2); slm.t

                %                 SITE_term = term(SITE_FN_disc);
                %                 M1 = 1 + SITE_term;
                %                 M2 = 1;
                %                 slm1  = SurfStatLinMod(bayes_validationvar_set_disc, M1);
                %                 slm2  = SurfStatLinMod(bayes_validationvar_set_disc, M2);
                %                 slm  = SurfStatF(slm1, slm2); slm.t
                
            end
            
        end
        
        zthickness_set_disc       = zscore(bayes_thickness_set_disc);
        zmyelin_set_disc          = zscore(bayes_myelin_set_disc);
        zclinicalvar_set_disc     = zscore(bayes_clinicalvar_set_disc);
        zvalidationvar_set_disc   = zscore(bayes_validationvar_set_disc);
        zclinicalvar2_set_disc     = zscore(bayes_clinicalvar2_set_disc);
        
    end
    
    % 2-2) replication process
    for replication_set = 1
        
        % 2-2-1) variable setup
        fieldnames_set        = LSPhenoDataSample_disc_repl.Properties.VariableNames;

        thickness_set_org_repl     = thickness_set_org(501:1000, :);
        myelin_set_org_repl        = myelin_set_org(501:1000, :);
        clinicalvar_set_org_repl   = phenotypic_data_replication(:, clinical_idx);
        clinicalvar_names_repl     = fieldnames_set(clinical_idx+1);
        clinicalvar2_set_org_repl  = cbcl_subscale_repl;
        clinicalvar2_names_repl    = cbclscore_name;
        validationvar_set_org_repl = phenotypic_data_replication(:, validation_idx);
        validationvar_names_repl   = fieldnames_set(validation_idx+1);
        demorace_data_org_repl     = demorace_data_org(501:1000, :);
        ethn_data_org_repl         = ethn_data_org(501:1000, :);
        parental_edu_org_repl      = parental_edu_org(501:1000, :);
        parental_income_org_repl   = parental_income_org(501:1000);
        houseincome_org_repl       = houseincome_org(501:1000);
        
        CASE_FN_repl       = CASE_org(501:1000);
        AGE_FN_repl        = AGE_org(501:1000);
        GENDER_FN_repl     = GENDER_org(501:1000);
        SITE_FN_repl       = cellstr(num2str(SITE_org(501:1000)));

        thickness_set_repl       = thickness_set_org_repl;
        myelin_set_repl          = myelin_set_org_repl;
        clinicalvar_set_repl     = clinicalvar_set_org_repl;
        validationvar_set_repl   = validationvar_set_org_repl;
        demorace_data_set_repl   = demorace_data_org_repl;
        ethn_data_set_repl       = ethn_data_org_repl;
        parental_edu_set_repl    = parental_edu_org_repl;
        parental_income_set_repl = parental_income_org_repl; 
        houseincome_set_repl     = houseincome_org_repl;
        clinicalvar2_set_repl    = clinicalvar2_set_org_repl;
        
        pheno_scores_repl = [ clinicalvar_set_repl validationvar_set_repl ];
        temp = corr(pheno_scores_repl);
        
        figure; heatmap(temp, fieldnames_set([clinical_idx validation_idx]+1), fieldnames_set([clinical_idx validation_idx]+1)); xtickangle(45); colormap('jet'); caxis([0 1]);
        
        % 2-1-2) outlier check
        threshold_for_outlier = 0.2;
        
        % - imaging data
        zthickness_set_org_repl = (thickness_set_org_repl - repmat(mean(thickness_set_org_repl, 1), size(thickness_set_org_repl, 1), 1))./repmat(std(thickness_set_org_repl, 0, 1), size(thickness_set_org_repl, 1), 1);
        thresholded_thickness_repl = abs(zthickness_set_org_repl)>3.09;
        thresholded_thickness_repl = sum(thresholded_thickness_repl, 2)/size(thresholded_thickness_repl, 2);
        [a b] = sort(thresholded_thickness_repl); a'
        outlier_idx_repl = unique(b(a>=threshold_for_outlier));
        
        zmyelin_set_org_repl = (myelin_set_org_repl - repmat(mean(myelin_set_org_repl, 1), size(myelin_set_org_repl, 1), 1))./repmat(std(myelin_set_org_repl, 0, 1), size(myelin_set_org_repl, 1), 1);
        thresholded_myelin_repl = abs(zmyelin_set_org_repl)>3.09;
        thresholded_myelin_repl = sum(thresholded_myelin_repl, 2)/size(thresholded_myelin_repl, 2);
        [a b] = sort(thresholded_myelin_repl); a'
        outlier_idx_repl = unique([ outlier_idx_repl; b(a>=threshold_for_outlier) ]);
        outlier_idx_repl = [ outlier_idx_repl; find(mean(zmyelin_set_org_repl(:, [ 12307 15046 1151 1945 ]), 2)>2)' ];
        
        thickness_set_repl(outlier_idx_repl, :)     = [];
        myelin_set_repl(outlier_idx_repl, :)        = [];
        clinicalvar_set_repl(outlier_idx_repl, :)   = [];
        validationvar_set_repl(outlier_idx_repl, :) = [];
        clinicalvar2_set_repl(outlier_idx_repl, :)  = [];
        CASE_FN_repl(outlier_idx_repl)              = [];
        AGE_FN_repl(outlier_idx_repl)               = [];
        GENDER_FN_repl(outlier_idx_repl)            = [];
        SITE_FN_repl(outlier_idx_repl)              = [];
        fn_prob_replication(outlier_idx_repl)      = [];
        demographic_data_replication(outlier_idx_repl, :) = [];
        demorace_data_set_repl(outlier_idx_repl, :) = [];
        ethn_data_set_repl(outlier_idx_repl) = [];
        parental_edu_set_repl(outlier_idx_repl)     = [];
        parental_income_set_repl(outlier_idx_repl)     = [];
        Release_replication(outlier_idx_repl)         = [];
        houseincome_set_repl(outlier_idx_repl)         = [];
        % 2-1-3) statistically remove the age, sex and site effects for brain imaging data and the site effect for phenotype scores
        for age_sex_correction = 1
            
            for thickness = 1
                
                AGE_term = term(AGE_FN_repl);
                SEX_term = term(GENDER_FN_repl);
                M = 1 + AGE_term + SEX_term;
                slm  = SurfStatLinMod(thickness_set_repl, M, ST);
                %                 slm  = SurfStatT(slm, AGE_FN_repl);
                %                 slm  = SurfStatT(slm, SEX_term.M-SEX_term.F);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(mean(thickness_set_repl, 1), ST);  SurfStatColLim([0 4.5]);
                
                thickness_set_repl_age_sex_corr = thickness_set_repl - slm.X(:, 2:4)*slm.coef(2:4, :);
                slm  = SurfStatLinMod(thickness_set_repl_age_sex_corr, M, ST);
                %                 slm  = SurfStatT(slm, AGE_FN_repl);
                %                 slm  = SurfStatT(slm, SEX_term.M-SEX_term.F);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(mean(thickness_set_repl_age_sex_corr, 1), ST);  SurfStatColLim([0 4.5]);
                
            end
            for myelin = 1
                
                AGE_term = term(AGE_FN_repl);
                SEX_term = term(GENDER_FN_repl);
                M = 1 + AGE_term + SEX_term;
                slm  = SurfStatLinMod(myelin_set_repl, M, ST);
                %                 slm  = SurfStatT(slm, AGE_FN_repl);
                %                 slm  = SurfStatT(slm, SEX_term.M-SEX_term.F);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(mean(myelin_set_repl, 1), ST);  SurfStatColLim([0.3 0.6]);
                
                myelin_set_repl_age_sex_corr = myelin_set_repl - slm.X(:, 2:4)*slm.coef(2:4, :);
                slm  = SurfStatLinMod(myelin_set_repl_age_sex_corr, M, ST);
                %                 slm  = SurfStatT(slm, AGE_FN_repl);
                %                 slm  = SurfStatT(slm, SEX_term.M-SEX_term.F);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(mean(myelin_set_repl_age_sex_corr, 1), ST); SurfStatColLim([0.25 0.35]);
                
            end
            
        end
        for site_correction = 1
            
            for thickness = 1
                
                thickness_set_std = std(thickness_set_repl_age_sex_corr, 0, 1);
                thickness_set_temp = thickness_set_repl_age_sex_corr;
                
                thickness_set_temp(:, find(thickness_set_std == 0)) = [];
                
                mod = [];
                bayes_thickness_set_temp = combat(thickness_set_temp', str2num(cell2mat(SITE_FN_repl)), mod, 1)';
                bayes_thickness_set_repl = zeros(length(CASE_FN_repl), 20484); bayes_thickness_set_repl(:,  find(thickness_set_std ~= 0)) = bayes_thickness_set_temp;
                
                SITE_term = term(SITE_FN_repl);
                M1 = 1 + SITE_term;
                M2 = 1;
                slm1  = SurfStatLinMod(thickness_set_repl_age_sex_corr, M1, ST);
                slm2  = SurfStatLinMod(thickness_set_repl_age_sex_corr, M2, ST);
                %                 slm  = SurfStatF(slm1, slm2);
                %                 [ pval, peak, clus, clusid ] = SurfStatP(slm, mask, 0.025);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(pval.C, ST); SurfStatColLim([0 0.05]);
                
                SITE_term = term(SITE_FN_repl);
                M1 = 1 + SITE_term;
                M2 = 1;
                slm1  = SurfStatLinMod(bayes_thickness_set_repl, M1, ST);
                slm2  = SurfStatLinMod(bayes_thickness_set_repl, M2, ST);
                %                 slm  = SurfStatF(slm1, slm2);
                %                 [ pval, peak, clus, clusid ] = SurfStatP(slm, mask, 0.025);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(pval.C, ST); SurfStatColLim([0 0.05]);
                
            end
            for myelin = 1
                
                myelin_set_std = std(myelin_set_repl_age_sex_corr, 0, 1);
                myelin_set_temp = myelin_set_repl_age_sex_corr;
                
                myelin_set_temp(:, find(myelin_set_std == 0)) = [];
                
                mod = [];
                bayes_myelin_set_temp = combat(myelin_set_temp', str2num(cell2mat(SITE_FN_repl)), mod, 1)';
                bayes_myelin_set_repl = zeros(length(CASE_FN_repl), 20484); bayes_myelin_set_repl(:,  find(myelin_set_std ~= 0)) = bayes_myelin_set_temp;
                
                SITE_term = term(SITE_FN_repl);
                M1 = 1 + SITE_term;
                M2 = 1;
                slm1  = SurfStatLinMod(myelin_set_repl_age_sex_corr, M1, ST);
                slm2  = SurfStatLinMod(myelin_set_repl_age_sex_corr, M2, ST);
                %                 slm  = SurfStatF(slm1, slm2);
                %                 [ pval, peak, clus, clusid ] = SurfStatP(slm, mask, 0.025);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(pval.C, ST); SurfStatColLim([0 0.05]);
                
                SITE_term = term(SITE_FN_repl);
                M1 = 1 + SITE_term;
                M2 = 1;
                slm1  = SurfStatLinMod(bayes_myelin_set_repl, M1, ST);
                slm2  = SurfStatLinMod(bayes_myelin_set_repl, M2, ST);
                %                 slm  = SurfStatF(slm1, slm2);
                %                 [ pval, peak, clus, clusid ] = SurfStatP(slm, mask, 0.025);
                %                 figure; SurfStatView(slm.t, ST); SurfStatColLim([-3 3]);
                %                 figure; SurfStatView(pval.C, ST); SurfStatColLim([0 0.05]);
                
            end
            for clinicalvar = 1
                
                clinicalvar_set_std = std(clinicalvar_set_repl, 0, 1);
                clinicalvar_set_temp = clinicalvar_set_repl;
                clinicalvar_set_temp(:, find(clinicalvar_set_std == 0)) = [];
                clinicalvar2_set_std = std(clinicalvar2_set_repl, 0, 1);
                clinicalvar2_set_temp = clinicalvar2_set_repl;
                clinicalvar2_set_temp(:, find(clinicalvar2_set_std == 0)) = [];
                
                mod = [];
                bayes_clinicalvar_set_temp = combat(clinicalvar_set_temp', str2num(cell2mat(SITE_FN_repl)), mod, 1)';
                bayes_clinicalvar_set_repl = zeros(length(CASE_FN_repl), size(clinicalvar_set_repl, 2)); bayes_clinicalvar_set_repl(:,  find(clinicalvar_set_std ~= 0)) = bayes_clinicalvar_set_temp;
                bayes_clinicalvar2_set_temp = combat(clinicalvar2_set_temp', str2num(cell2mat(SITE_FN_repl)), mod, 1)';
                bayes_clinicalvar2_set_repl = zeros(length(CASE_FN_repl), size(clinicalvar2_set_repl, 2)); bayes_clinicalvar2_set_repl(:,  find(clinicalvar2_set_std ~= 0)) = bayes_clinicalvar2_set_temp;
                %                 SITE_term = term(SITE_FN_repl);
                %                 M1 = 1 + SITE_term;
                %                 M2 = 1;
                %                 slm1  = SurfStatLinMod(clinicalvar_set_repl, M1);
                %                 slm2  = SurfStatLinMod(clinicalvar_set_repl, M2);
                %                 slm  = SurfStatF(slm1, slm2); slm.t

                %                 SITE_term = term(SITE_FN_repl);
                %                 M1 = 1 + SITE_term;
                %                 M2 = 1;
                %                 slm1  = SurfStatLinMod(bayes_clinicalvar_set_repl, M1);
                %                 slm2  = SurfStatLinMod(bayes_clinicalvar_set_repl, M2);
                %                 slm  = SurfStatF(slm1, slm2); slm.t
                
            end
            for validationlvar = 1
                
                validationvar_set_std = std(validationvar_set_repl, 0, 1);
                validationvar_set_temp = validationvar_set_repl;
                validationvar_set_temp(:, find(validationvar_set_std == 0)) = [];
                
                mod = [];
                bayes_validationvar_set_temp = combat(validationvar_set_temp', str2num(cell2mat(SITE_FN_repl)), mod, 1)';
                bayes_validationvar_set_repl = zeros(length(str2num(cell2mat(SITE_FN_repl))), size(validationvar_set_repl, 2)); bayes_validationvar_set_repl(:,  find(validationvar_set_std ~= 0)) = bayes_validationvar_set_temp;
                %                 SITE_term = term(SITE_FN_repl);
                %                 M1 = 1 + SITE_term;
                %                 M2 = 1;
                %                 slm1  = SurfStatLinMod(validationvar_set_repl, M1);
                %                 slm2  = SurfStatLinMod(validationvar_set_repl, M2);
                %                 slm  = SurfStatF(slm1, slm2); slm.t

                %                 SITE_term = term(SITE_FN_repl);
                %                 M1 = 1 + SITE_term;
                %                 M2 = 1;
                %                 slm1  = SurfStatLinMod(bayes_validationvar_set_repl, M1);
                %                 slm2  = SurfStatLinMod(bayes_validationvar_set_repl, M2);
                %                 slm  = SurfStatF(slm1, slm2); slm.t
                
            end
            
        end
        
        zthickness_set_repl       = zscore(bayes_thickness_set_repl);
        zmyelin_set_repl          = zscore(bayes_myelin_set_repl);
        zclinicalvar_set_repl     = zscore(bayes_clinicalvar_set_repl);
        zvalidationvar_set_repl   = zscore(bayes_validationvar_set_repl);
        zclinicalvar2_set_repl     = zscore(bayes_clinicalvar2_set_repl);
        
    end
    
    % 2-3) demographic comparison
    for demo_comp = 1
        
        % age
        [ mean(AGE_FN_disc) std(AGE_FN_disc);
            mean(AGE_FN_repl) std(AGE_FN_repl); ]
        
        [ h p ci stats ] = ttest2(AGE_FN_disc, AGE_FN_repl)
        
        % sex
        temp = [ sum(strcmp(GENDER_FN_disc, 'M')) sum(strcmp(GENDER_FN_disc, 'F'));
                 sum(strcmp(GENDER_FN_repl, 'M')) sum(strcmp(GENDER_FN_repl, 'F')) ]
     
        [h,p,X2] = chi2cont(temp)
        
        % site
        temp1 = tabulate(SITE_FN_disc);
        temp2 = tabulate(SITE_FN_repl);
        [a b] = sort(temp1(:, 1));
        [c d] = sort(temp2(:, 1));
        temp = [ cell2mat(temp1(b, 2))'; cell2mat(temp2(d, 2))' ];
        [h,p,X2] = chi2cont(temp)
        
        % race
        demorace_data_set_disc_unirace = demorace_data_set_disc;
        demorace_data_set_disc_unirace(find(sum(demorace_data_set_disc, 2)>1), :) = [];
        demorace_data_set_disc_multirace = demorace_data_set_disc(find(sum(demorace_data_set_disc, 2)>1), :);
        demorace_data_set_repl_unirace = demorace_data_set_repl;
        demorace_data_set_repl_unirace(find(sum(demorace_data_set_repl, 2)>1), :) = [];
        demorace_data_set_repl_multirace = demorace_data_set_repl(find(sum(demorace_data_set_repl, 2)>1), :);
        
        temp_unirace = [ sum(demorace_data_set_disc_unirace); sum(demorace_data_set_repl_unirace) ]
        [h,p,X2] = chi2cont(temp_unirace)
        
        W_total   = temp_unirace(:, 1);
        B_total   = temp_unirace(:, 2);
        A_total   = sum(temp_unirace(:, 9:15), 2);
        O_total   = sum(temp_unirace(:, [ 3:8 16 ]), 2);
        NR_total  = sum(temp_unirace(:, [17 18]), 2);
        MR_total  = [ size(demorace_data_set_disc_multirace, 1); size(demorace_data_set_repl_multirace, 1) ];
        
        temp2 = [ W_total B_total  A_total O_total MR_total NR_total ];
        [h,p,X2] = chi2cont(temp2)
        
        % ethnicity
        temp1 = tabulate(ethn_data_set_disc);
        temp2 = tabulate(ethn_data_set_repl);
        [ temp1(temp1(:, 2) ~= 0, 2) temp2(temp2(:, 2) ~= 0, 2) ]
        p=myfisher([ temp1(temp1(:, 2) ~= 0, 2) temp2(temp2(:, 2) ~= 0, 2) ])
        
        % education level of parents
        temp = tabulate(parental_edu_set_disc);
        temp1 = sum(temp(temp(:, 1)<15, 2))                 % E/M/H
        temp2 = sum(temp(temp(:, 1)>=15&temp(:, 1)<19, 2))  % B
        temp3 = sum(temp(temp(:, 1)>=19, 2))                % G
        
        temp = tabulate(parental_edu_set_repl);
        temp4 = sum(temp(temp(:, 1)<15, 2))                 % E/M/H
        temp5 = sum(temp(temp(:, 1)>=15&temp(:, 1)<19, 2))  % B
        temp6 = sum(temp(temp(:, 1)>=19, 2))                % G
        [h,p,X2] = chi2cont([ temp1 temp2 temp3; temp4 temp5 temp6 ])
        
        % House income of parents
        temp1 = sum(houseincome_set_disc<=6)
        temp2 = sum(houseincome_set_disc>6&houseincome_set_disc<=8)
        temp3 = sum(houseincome_set_disc>8&houseincome_set_disc<11)
        
        temp4 = sum(houseincome_set_repl<=6)
        temp5 = sum(houseincome_set_repl>6&houseincome_set_repl<=8)
        temp6 = sum(houseincome_set_repl>8&houseincome_set_repl<11)
        [h,p,X2] = chi2cont([ temp1 temp2 temp3; temp4 temp5 temp6 ])
        
    end
    
    
    FPATH = '/data3/cdb/shong/01_project/02_abcd_brain_development/03_result/';
    %     T = array2table(zthickness_set_disc, 'RowNames', CASE_FN_disc);
    %     writetable(T, 'CMI_ABCD_Sample1_thickness_disc.csv', 'WriteRowNames', true);
    %     T = array2table(zmyelin_set_disc, 'RowNames', CASE_FN_disc);
    %     writetable(T, 'CMI_ABCD_Sample1_myelin_disc.csv',  'WriteRowNames', true);
    %     T = array2table(zclinicalvar_set_disc, 'RowNames', CASE_FN, 'VariableNames', clinicalvar_names_disc);
    %     writetable(T, 'CMI_ABCD_Sample1_clinicalvars_disc.csv',  'WriteRowNames', true);
    %     T = array2table(zvalidationvar_set_disc, 'RowNames', CASE_FN, 'VariableNames', validationvar_names_disc);
    %     writetable(T, 'CMI_ABCD_Sample1_valdiationvar_disc.csv',  'WriteRowNames', true);
    
    %     T = array2table(zthickness_set_repl, 'RowNames', CASE_FN_repl);
    %     writetable(T, 'CMI_ABCD_Sample1_thickness_repl.csv', 'WriteRowNames', true);
    %     T = array2table(zmyelin_set_repl, 'RowNames', CASE_FN_repl);
    %     writetable(T, 'CMI_ABCD_Sample1_myelin_repl.csv',  'WriteRowNames', true);
    %     T = array2table(zclinicalvar_set_repl, 'RowNames', CASE_FN, 'VariableNames', clinicalvar_names_repl);
    %     writetable(T, 'CMI_ABCD_Sample1_clinicalvars_repl.csv',  'WriteRowNames', true);
    %     T = array2table(zvalidationvar_set_repl, 'RowNames', CASE_FN, 'VariableNames', validationvar_names_repl);
    %     writetable(T, 'CMI_ABCD_Sample1_valdiationvar_repl.csv',  'WriteRowNames', true);
    
    %     figure; SurfStatView(mean(bayes_thickness_set_disc, 1), ST); SurfStatColLim([0 3]);
    %     export_fig([FPATH '/99_thickness_mean_discovery' ], '-m2', '-png');  close(gcf);
    %     figure; SurfStatView(mean(bayes_myelin_set_disc, 1), ST); SurfStatColLim([0.2 0.35]); colormap('gray');
    %     export_fig([FPATH '/99_myelin_mean_discovery' ], '-m2', '-png');  close(gcf);
    %     figure; SurfStatView(mean(bayes_thickness_set_repl, 1), ST); SurfStatColLim([0 3]);
    %     export_fig([FPATH '/99_thickness_mean_replication' ], '-m2', '-png');  close(gcf);
    %     figure; SurfStatView(mean(bayes_myelin_set_repl, 1), ST); SurfStatColLim([0.2 0.35]); colormap('gray');
    %     export_fig([FPATH '/99_myelin_mean_replication' ], '-m2', '-png');  close(gcf);
    
end

%% 3) Similarity Network Fusion (SNF) analysis
for SNFanalysis = 1
    
    for imaging_network = 1
        
       imaging_features_disc = [ zthickness_set_disc(:, ~any(isnan(zthickness_set_disc))) zmyelin_set_disc(:, ~any(isnan(zmyelin_set_disc))) ]; 
       imaging_features_repl = [ zthickness_set_repl(:, ~any(isnan(zthickness_set_repl))) zmyelin_set_repl(:, ~any(isnan(zmyelin_set_repl))) ];      
 
    end
    
    for clinicalvar_network = 1
        
        clinical_features_disc = zclinicalvar_set_disc(:,  ~any(isnan(zclinicalvar_set_disc)));
        clinical_features_repl = zclinicalvar_set_repl(:,  ~any(isnan(zclinicalvar_set_repl)));
        
        validate_features_disc = zvalidationvar_set_disc(:,  ~any(isnan(zvalidationvar_set_disc)));
        validate_features_repl = zvalidationvar_set_repl(:,  ~any(isnan(zvalidationvar_set_repl)));
        
    end
    
    for SNF_analysis = 1
        
        %%% default parameter setting
        K     = 30;  %number of neighbors, usually (10~30)
        alpha = 0.5; %hyperparameter, usually (0.3~0.8)
        T     = 20;  %Number of Iterations, usually (10~30)
        iteration_num = 100;
        rsamplerate = 0.9;
        candidate_K_num = 20;
        warning off;
        
        howmanytimesselected_disc = zeros(length(fn_prob_discovery), length(fn_prob_discovery));
        howmanytimesselected_repl = zeros(length(fn_prob_replication), length(fn_prob_replication));
        consensus_matrix_disc = zeros(length(fn_prob_discovery), length(fn_prob_discovery), candidate_K_num);
        consensus_matrix_repl = zeros(length(fn_prob_replication), length(fn_prob_replication), candidate_K_num);
        for i = 1 : iteration_num
            
            i
            y_disc = sort(datasample(1:length(fn_prob_discovery), round(length(fn_prob_discovery)*rsamplerate), 'Replace', false)); %, 'Weights', fn_prob_discovery);
            y_repl = sort(datasample(1:length(fn_prob_replication), round(length(fn_prob_replication)*rsamplerate), 'Replace', false)); %, 'Weights', fn_prob_replication);
            
            tempdisc = nchoosek(y_disc, 2);
            tempind1 = sub2ind([length(fn_prob_discovery), length(fn_prob_discovery)], tempdisc(:, 1), tempdisc(:, 2));
            tempind2 = sub2ind([length(fn_prob_discovery), length(fn_prob_discovery)], tempdisc(:, 2), tempdisc(:, 1));
            howmanytimesselected_disc(tempind1) = howmanytimesselected_disc(tempind1) + 1; 
            howmanytimesselected_disc(tempind2) = howmanytimesselected_disc(tempind2) + 1; 
            temprepl = nchoosek(y_repl, 2);
            tempind1 = sub2ind([length(fn_prob_replication), length(fn_prob_replication)], temprepl(:, 1), temprepl(:, 2));
            tempind2 = sub2ind([length(fn_prob_replication), length(fn_prob_replication)], temprepl(:, 2), temprepl(:, 1));
            howmanytimesselected_repl(tempind1) = howmanytimesselected_repl(tempind1) + 1; 
            howmanytimesselected_repl(tempind2) = howmanytimesselected_repl(tempind2) + 1; 
            
%             [coeff_disc, score_disc, latent_disc, tsquared_disc, explained_disc, mu_disc] = pca(imaging_features_disc(y_disc, :));
%             [coeff_repl, score_repl, latent_repl, tsquared_repl, explained_repl, mu_repl] = pca(imaging_features_repl(y_repl, :));
%             
%             compnum_disc = min(find(cumsum(explained_disc)>95));
%             compnum_repl = min(find(cumsum(explained_repl)>95));
            
%             figure; imagesc(score_disc(:, 1:compnum_disc)); colormap('jet'); colorbar; caxis([-10 10]);
%             figure; imagesc(score_repl(:, 1:compnum_repl)); colormap('jet'); colorbar; caxis([-10 10]);
            score_disc = imaging_features_disc(y_disc, :); compnum_disc = 40968;
            score_repl = imaging_features_repl(y_repl, :); compnum_repl = 40968;
            
            Data1_disc = Standard_Normalization(score_disc(:, 1:compnum_disc));
            Data2_disc = Standard_Normalization(clinical_features_disc(y_disc, :));
            Data1_repl = Standard_Normalization(score_repl(:, 1:compnum_repl));
            Data2_repl = Standard_Normalization(clinical_features_repl(y_repl, :));
          
            %%%Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows;
            Data1_disc = squareform(pdist(Data1_disc, 'cosine'));
            Data2_disc = squareform(pdist(Data2_disc, 'cosine'));
            Data1_repl = squareform(pdist(Data1_repl, 'cosine'));
            Data2_repl = squareform(pdist(Data2_repl, 'cosine'));
            
            %%%next, construct similarity graphs
            W1_disc = affinityMatrix(Data1_disc, K, alpha);
            W2_disc = affinityMatrix(Data2_disc, K, alpha);
            W1_repl = affinityMatrix(Data1_repl, K, alpha);
            W2_repl = affinityMatrix(Data2_repl, K, alpha);
            
            %next, we fuse all the graphs
            % then the overall matrix can be computed by similarity network fusion(SNF):
            W_disc = SNF({W1_disc,W2_disc}, K, T);
            W_repl = SNF({W1_repl,W2_repl}, K, T);
            
            %%%%With this unified graph W of size n x n, you can do either spectral clustering or Kernel NMF. If you need help with further clustering, please let us know.
            %%for example, spectral clustering
            for C = 2 : candidate_K_num
            
                group_disc = SpectralClustering(W_disc,C); %%%the final subtypes information
                group_repl = SpectralClustering(W_repl,C); %%%the final subtypes information
                
                consensus_matrix_disc_temp = zeros(length(fn_prob_discovery), length(fn_prob_discovery));
                consensus_matrix_repl_temp = zeros(length(fn_prob_replication), length(fn_prob_replication));
                for j = 1 : C
                   
                    % discovery
                    if(length(find(group_disc == j)) == 1)
                        comb = [ find(group_disc == j) find(group_disc == j) ];
                    else
                        comb = nchoosek(find(group_disc == j), 2);
                        comb = [ comb; find(group_disc == j) find(group_disc == j) ];
                    end
                    indices1 = sub2ind([length(fn_prob_discovery), length(fn_prob_discovery)], y_disc(comb(:, 1))', y_disc(comb(:, 2))');
                    indices2 = sub2ind([length(fn_prob_discovery), length(fn_prob_discovery)], y_disc(comb(:, 2))', y_disc(comb(:, 1))');
                    consensus_matrix_disc_temp(indices1) = 1;
                    consensus_matrix_disc_temp(indices2) = 1;
                    
                    % replication
                    if(length(find(group_repl == j)) == 1)
                        comb = [ find(group_repl == j) find(group_repl == j) ];
                    else
                        comb = nchoosek(find(group_repl == j), 2);
                        comb = [ comb; find(group_repl == j) find(group_repl == j) ];
                    end
                    indices1 = sub2ind([length(fn_prob_replication), length(fn_prob_replication)], y_repl(comb(:, 1))', y_repl(comb(:, 2))');
                    indices2 = sub2ind([length(fn_prob_replication), length(fn_prob_replication)], y_repl(comb(:, 2))', y_repl(comb(:, 1))');
                    consensus_matrix_repl_temp(indices1) = 1;
                    consensus_matrix_repl_temp(indices2) = 1;
                    
                end
                consensus_matrix_disc(:, :, C) = consensus_matrix_disc(:, :, C) + consensus_matrix_disc_temp;
                consensus_matrix_repl(:, :, C) = consensus_matrix_repl(:, :, C) + consensus_matrix_repl_temp;
                
            end
            
        end
        
        Ci_disc = [];
        Q_disc  = [];
        Ci_repl = [];
        Q_repl  = [];
        maxcomm_disc = [];
        maxcomm_repl = [];
        for C = 2 : candidate_K_num
            
            temp1 = consensus_matrix_disc(:, :, C)./howmanytimesselected_disc; temp1(logical(eye(size(temp1)))) = 1;
            temp2 = consensus_matrix_repl(:, :, C)./howmanytimesselected_repl; temp2(logical(eye(size(temp2)))) = 1;
            [ group_disc, ~ ] = SpectralClustering(temp1, C);
            [ group_repl, ~ ] = SpectralClustering(temp2, C);
            [ Ci_disc(:, C-1) Q_disc(C-1) ] = community_louvain(temp1, [], group_disc);
            [ Ci_repl(:, C-1) Q_repl(C-1) ] = community_louvain(temp2, [], group_repl);
            maxcomm_disc(C-1) = max(Ci_disc(:, C-1));
            maxcomm_repl(C-1) = max(Ci_repl(:, C-1));
            
        end
        
        figure; plot(Q_disc);
        figure; plot(Q_repl);
        curridx = 6; [On,Wr] = reorder_mod(consensus_matrix_disc(:, :, curridx), Ci_disc(:, curridx-1)); 
        figure; imagesc(Wr); colorbar; colormap('jet'); caxis([0 iteration_num]);
        curridx = 6; [On,Wr] = reorder_mod(consensus_matrix_repl(:, :, curridx), Ci_repl(:, curridx-1)); 
        figure; imagesc(Wr); colorbar; colormap('jet'); caxis([0 iteration_num]);
        
        for Monti2003_consensus_statistics = 1
            
            legend_str = cell(1, candidate_K_num - 1);
            for C = 2 : candidate_K_num
                legend_str{C-1} = num2str(C);
            end
            
            A_disc = zeros(1, candidate_K_num - 2 + 1);
            A_repl = zeros(1, candidate_K_num - 2 + 1);
            CDF_consensus_disc = zeros(candidate_K_num - 2 + 1, length(0 : 0.01 : 1));
            CDF_consensus_repl = zeros(candidate_K_num - 2 + 1, length(0 : 0.01 : 1));
            for C = 2 : candidate_K_num
                
                temp1 = consensus_matrix_disc(:, :, C)./howmanytimesselected_disc; temp1(logical(eye(size(temp1)))) = 0;
                temp1 = tril(temp1); temp1 = temp1(temp1~=0);
                
                count = 1;
                for k = 0 : 0.01 : 1
                    CDF_consensus_disc(C-1, count) = sum(temp1 <= k)/(size(consensus_matrix_disc, 1)*(size(consensus_matrix_disc, 1)-1)/2);
                    count = count + 1;
                end
                A_disc(C-1) = sum(CDF_consensus_disc(C-1, :)*0.01);
                
                temp2 = consensus_matrix_repl(:, :, C)./howmanytimesselected_repl; temp2(logical(eye(size(temp2)))) = 1;
                temp2 = tril(temp2); temp2 = temp2(temp2~=0);
                
                count = 1;
                for k = 0 : 0.01 : 1
                    CDF_consensus_repl(C-1, count) = sum(temp2 <= k)/(size(consensus_matrix_repl, 1)*(size(consensus_matrix_repl, 1)-1)/2);
                    count = count + 1;
                end
                A_repl(C-1) = sum(CDF_consensus_repl(C-1, :)*0.01);
                
            end
            
            delA_disc = zeros(1, candidate_K_num);
            delA_repl = zeros(1, candidate_K_num);
            for C = 2 : candidate_K_num-1
                
                if(C==2)
                    delA_disc(C) = A_disc(C-1);
                else
                    delA_disc(C) = (A_disc(C-1+1) -  A_disc(C-1))/A_disc(C-1);
                end
                if(C==2)
                    delA_repl(C) = A_repl(C-1);
                else
                    delA_repl(C) = (A_repl(C-1+1) -  A_repl(C-1))/A_repl(C-1);
                end
                
            end
            
            figure; plot(delA_disc(2:candidate_K_num));
            figure; plot(delA_repl(2:candidate_K_num));
            
        end
        
        Ci_disc = [];
        Ci_repl = [];
        within_conavg_disc = [];
        between_conavg_disc = [];
        within_conavg_repl = [];
        between_conavg_repl = [];
        for C = 2 : candidate_K_num
            
            temp1 = consensus_matrix_disc(:, :, C)./howmanytimesselected_disc; temp1(logical(eye(size(temp1)))) = 1;
            temp2 = consensus_matrix_repl(:, :, C)./howmanytimesselected_repl; temp2(logical(eye(size(temp2)))) = 1;
            [ Ci_disc(:, C-1), ~ ] = SpectralClustering(temp1, C);
            [ Ci_repl(:, C-1), ~ ] = SpectralClustering(temp2, C);
            
            within_comb = [];
            within_consensus = [];
            for k = 1 : C
               
                within_idx = find(Ci_disc(:, C-1) == k);
                comb = nchoosek(within_idx, 2);
                comb1 = sub2ind([ size(temp1, 1) size(temp1, 2) ], comb(:, 1), comb(:, 2));
                comb2 = sub2ind([ size(temp1, 1) size(temp1, 2) ], comb(:, 2), comb(:, 1));
                within_comb = [ within_comb comb1' comb2' ];
                within_consensus = [ within_consensus temp1([ comb1' comb2' ]) ]; 
                
            end
            
            total_comb = 1 : size(temp1, 1)*size(temp1, 2);
            total_comb(sort(within_comb)) = [];
            between_consensus = temp1(total_comb);
            within_conavg_disc = [ within_conavg_disc mean(within_consensus) ];
            between_conavg_disc = [ between_conavg_disc mean(between_consensus) ];
            
            within_comb = [];
            within_consensus = [];
            for k = 1 : C
               
                within_idx = find(Ci_repl(:, C-1) == k);
                comb = nchoosek(within_idx, 2);
                comb1 = sub2ind([ size(temp2, 1) size(temp2, 2) ], comb(:, 1), comb(:, 2));
                comb2 = sub2ind([ size(temp2, 1) size(temp2, 2) ], comb(:, 2), comb(:, 1));
                within_comb = [ within_comb comb1' comb2' ];
                within_consensus = [ within_consensus temp2([ comb1' comb2' ]) ]; 
                
            end
            
            total_comb = 1 : size(temp2, 1)*size(temp2, 2);
            total_comb(sort(within_comb)) = [];
            between_consensus = temp2(total_comb);
            within_conavg_repl = [ within_conavg_repl mean(within_consensus) ];
            between_conavg_repl = [ between_conavg_repl mean(between_consensus) ];
            
        end
        
        consensus_statistics_disc = within_conavg_disc./between_conavg_disc./(2:20);
        consensus_statistics_repl = within_conavg_repl./between_conavg_repl./(2:20);
        
        figure; plot(consensus_statistics_disc);
        figure; plot(consensus_statistics_repl);
        
        curridx = 2; 
        
        figure; 
        [On,Wr] = reorder_mod(consensus_matrix_disc(:, :, curridx), Ci_disc(:, curridx-1)); 
        lineidx = find(diff(Ci_disc(On, curridx-1)));
        subplot(1, 2, 1); imagesc(Wr); colorbar; colormap('jet'); caxis([0 iteration_num]);
        for i = 1 : length(lineidx)
            
            l = line([lineidx(i) lineidx(i)], [0, size(consensus_matrix_disc, 1)]); l.LineWidth = 3; l.Color = 'w';
            l = line([0, size(consensus_matrix_disc, 1)], [lineidx(i) lineidx(i)]); l.LineWidth = 3; l.Color = 'w';
            
        end
        [On,Wr] = reorder_mod(consensus_matrix_repl(:, :, curridx), Ci_repl(:, curridx-1));
        lineidx = find(diff(Ci_repl(On, curridx-1)));
        subplot(1, 2, 2); imagesc(Wr); colorbar; colormap('jet'); caxis([0 iteration_num]);
        for i = 1 : length(lineidx)
            
            l = line([lineidx(i) lineidx(i)], [0, size(consensus_matrix_disc, 1)]); l.LineWidth = 3; l.Color = 'w';
            l = line([0, size(consensus_matrix_disc, 1)], [lineidx(i) lineidx(i)]); l.LineWidth = 3; l.Color = 'w';
            
        end
        set(gcf, 'Position', [ -1411 234 1269 426 ]);
        
        tabulate(Ci_disc(:, curridx-1));
        tabulate(Ci_repl(:, curridx-1));
        
        curridx = 3;
        
        % discovery
        subtypes = Ci_disc(:, curridx-1);
        tabulate(subtypes)
        thickness_pattern = zeros(curridx, 20484);
        myelin_pattern    = zeros(curridx, 20484);
        for i = 1 : curridx
        
            thickness_pattern(i, :) = mean(imaging_features_disc(subtypes==i, 1:20484), 1);
            myelin_pattern(i, :)    = mean(imaging_features_disc(subtypes==i, 20485:40968), 1);
            
%             thickness_pattern(i, :) = mean(bayes_thickness_set_disc(subtypes==i, :), 1);
%             myelin_pattern(i, :)    = mean(bayes_myelin_set_disc(subtypes==i, :), 1);
            
        end
        figure; CSFSurfStatView(thickness_pattern, ST, 'INC', [ num2str(curridx) 'Dsubtypes thickness' ]); % CSFSurfStatViewColLim([1 2.6]);
        figure; CSFSurfStatView(myelin_pattern, ST, 'INC', [ num2str(curridx) 'Dsubtypes myelin' ]);  % CSFSurfStatViewColLim([0.25 0.35]);
        
        % replication
        subtypes = Ci_repl(:, curridx-1);
        tabulate(subtypes)
        thickness_pattern = zeros(curridx, 20484);
        myelin_pattern    = zeros(curridx, 20484);
        for i = 1 : curridx
        
            thickness_pattern(i, :) = mean(imaging_features_repl(subtypes==i, 1:20484), 1);
            myelin_pattern(i, :)    = mean(imaging_features_repl(subtypes==i, 20485:40968), 1);
            
%             thickness_pattern(i, :) = mean(bayes_thickness_set_repl(subtypes==i, :), 1);
%             myelin_pattern(i, :)    = mean(bayes_myelin_set_repl(subtypes==i, :), 1);
            
        end
        figure; CSFSurfStatView(thickness_pattern, ST, 'INC', [ num2str(curridx) 'Rsubtypes thickness' ]); % CSFSurfStatViewColLim([1 2.6]);
        figure; CSFSurfStatView(myelin_pattern, ST, 'INC', [ num2str(curridx) 'Rsubtypes myelin' ]); % CSFSurfStatViewColLim([0.25 0.35]);
         
        for i = 1 : iteration_num
           
            indices_disc = crossvalind('Kfold',1:length(fn_prob_discovery), 10);
            
            for C = 2 : candidate_K_num
                
                
                
            end
            
        end
        
    end
    
    for subtype_profiling = 1
        
        % 1 = Strongly Disagree; 2 = Disagree; 3 = Neutral (neither agree nor disagree); 4 = Agree; 5 = Strongly Agree
        
        % Environmental factors
        % meim_p_ss_exp: MEIM-R (Multigroup Ethnic Identity Measure – Revised) Exploration Subscale, Mean: (meim_1_p + meim_4_p + meim_5_p)/3;
        %       meim_1_p=I have spent time trying to find out more about my ethnic group, such as its history, traditions, and customs.
        %       meim_4_p=I have often done things that will help me understand my ethnic background better.
        %       meim_5_p=I have often talked to other people in order to learn more about my ethnic group.
        
        % meim_p_ss_com: MEIM-R Commitment and Attachment Subscale, Mean: (meim_2_p + meim_3_p + meim_6_p)/3;
        %      meim_2_p=I have a strong sense of belonging to my own ethnic group.
        %      meim_3_p=I understand pretty well what my ethnic group membership means to me.
        %      meim_6_p=I feel a strong attachment towards my own ethnic group.
        
        % meim_p_ss_total: MEIM-R Overall Mean: (meim_1_p + meim_2_p + meim_3_p + meim_4_p + meim_5_p + meim_6_p)/6
        
        % nsc_p_ss_mean_3_items: Neighborhood Safety Protocol, Mean of Parent Report, (neighborhood1r_p + neighborhood2r_p + neighborhood3r_p)/3;
        %       neighborhood1r_p=I feel safe walking in my neighborhood, day or night.
        %       neighborhood2r_p=Violence is not a problem in my neighborhood.
        %       neighborhood3r_p=My neighborhood is safe from crime.
        
        % fes_p_ss_fc: Conflict subscale from the Family Environment Scale Sum of Parent Report: fam_enviro1_p + fam_enviro2r_p + fam_enviro3_p + fam_enviro4r_p + fam_enviro5_p + fam_enviro6_p + fam_enviro7r_p + fam_enviro8_p + fam_enviro9r_p;
        %       fam_enviro1_p=We fight a lot in our family.
        %       fam_enviro2_p=Family members rarely become openly angry.
        %       fam_enviro3_p=Family members sometimes get so angry they throw things.
        %       fam_enviro4_p=Family members hardly ever lose their tempers.
        %       fam_enviro5_p=Family members often criticize each other.
        %       fam_enviro6_p=Family members sometimes hit each other.
        %       fam_enviro7_p=If there is a disagreement in our family, we try hard to smooth things over and keep the peace.
        %       fam_enviro8_p=Family members often try to one-up or outdo each other.
        %       fam_enviro9_p=In our family, we believe you don't ever get anywhere by raising your voice.
        
        % fes_y_ss_fc: Conflict Subscale from the Family Environment Scale Sum of Youth Report : fes_youth_q1 + fes_youth_q2 + fes_youth_q3 + fes_youth_q4 + fes_youth_q5 + fes_youth_q6 + fes_youth_q7 + fes_youth_q8 + fes_youth_q9
        %       fes_youth_q*=Same for parent version
        
        % pmq_y_ss_mean: Parental Monitoring: Mean (parent_monitor_q1_y + parent_monitor_q2_y + parent_monitor_q3_y + parent_monitor_q4_y parent_monitor_q5_y)/5
        %      parent_monitor_q1_y=How often do your parents/guardians know where you are?
        %      parent_monitor_q2_y=How often do your parents know who you are with when you are not at school and away from home?
        %      parent_monitor_q3_y=If you are at home when your parents or guardians are not, how often do you know how to get in touch with them?
        %      parent_monitor_q4_y=How often do you talk to your mom/dad or guardian about your plans for the coming day, such as your plans about what will happen at school or what you are going to do with friends?
        %      parent_monitor_q5_y=In an average week, how many times do you and your parents/guardians, eat dinner together?
        
        % crpbi_y_ss_parent: Acceptance Subscale Mean of Report by Parent Completing Protocol by youth: (crpbi_parent1_y + crpbi_parent2_y + crpbi_parent3_y], crpbi_parent4_y + crpbi_parent5_y)/5
        %      crpbi_parent1_y=First caregiver (caregiver participating in study/completing protocol). Makes me feel better after talking over my worries with him/her
        %      crpbi_parent2_y=First caregiver (caregiver participating in study/completing protocol). Smiles at me very often.
        %      crpbi_parent3_y=First caregiver (caregiver participating in study/completing protocol). Is able to make me feel better when I am upset.
        %      crpbi_parent4_y=First caregiver (caregiver participating in study/completing protocol). Believes in showing his/her love for me.
        %      crpbi_parent5_y=First caregiver (caregiver participating in study/completing protocol). Is easy to talk to.
        
        % crpbi_y_ss_caregiver: CRPBI - Acceptance Subscale Mean of Report by Secondary Caregiver by youth: (crpbi_caregiver12_y + crpbi_caregiver13_y + crpbi_caregiver14_y + crpbi_caregiver15_y + crpbi_caregiver16_y)/5
        %      crpbi_caregiver12_y=Second caregiver. Makes me feel better after talking over my worries with him/her.
        %      crpbi_caregiver13_y=Second caregiver. Smiles at me very often.
        %      crpbi_caregiver14_y=Second caregiver. Is able to make me feel better when I am upset.
        %      crpbi_caregiver15_y=Second caregiver. Believes in showing his/her love for me.
        %      crpbi_caregiver16_y=Second caregiver. Is easy to talk to.
        
        % srpf_y_ss_ses: SRPF School Environment Subscale, Sum: school_2_y + school_3_y + school_4_y + school_5_y + school_6_y + school_7_y
        %      school_2_y=In my school, students have lots of chances to help decide things like class activities and rules.
        %      school_3_y=I get along with my teachers.
        %      school_4_y=My teacher(s) notices when I am doing a good job and lets me know about it.
        %      school_5_y=There are lots of chances for students in my school to get involved in sports, clubs, or other school activities outside of class.
        %      school_6_y=I feel safe at my school.
        %      school_7_y=The school lets my parents know when I have done something well.
        
        % srpf_y_ss_iiss: SRPF School Involvement Subscale, Sum: school_8_y + school_9_y + school_10_y + school_12_y;
        %      school_8_y=I like school because I do well in class.
        %      school_9_y=I feel I'm just as smart as other kids my age.
        %      school_10_y=There are lots of chances to be part of class discussions or activities.
        %      school_12_y=In general, I like school a lot.
        
        % srpf_y_ss_dfs: SRPF School Disengagement Subscale, Sum: school_15_y + school_17_y;
        %      school_15_y=Usually, school bores me.
        %      school_17_y=Getting good grades is not so important to me.
        
        % trauma_num
        
        % Validation
        % cbcl_scr_syn_internal_t
        % cbcl_scr_syn_external_t
        % cbcl_scr_syn_totprob_t
        
        varNames = {'age', 'sex'}; 
        varNames((end+1):(end+length(clinicalvar_names))) = clinicalvar_names; 
        varNames((end+1):(end+length(validationvar_names))) = validationvar_names; 
        varNames(end+1) = {'subtype'};
        
        democlinical_table = table(AGE_FN/12, strcmp(GENDER_FN, 'F'), ...
            bayes_clinicalvar_set(:, 1), ...
            bayes_clinicalvar_set(:, 2), ...
            bayes_clinicalvar_set(:, 3), ...
            bayes_clinicalvar_set(:, 4), ...
            bayes_clinicalvar_set(:, 5), ...
            bayes_clinicalvar_set(:, 6), ...
            bayes_clinicalvar_set(:, 7), ...
            bayes_clinicalvar_set(:, 8), ...
            bayes_clinicalvar_set(:, 9), ...
            bayes_clinicalvar_set(:, 10), ...
            bayes_clinicalvar_set(:, 11), ...
            bayes_clinicalvar_set(:, 12), ...
            bayes_clinicalvar_set(:, 13), ...
            bayes_validationvar_set(:, 1), ...
            bayes_validationvar_set(:, 2), ...
            bayes_validationvar_set(:, 3), ...
            group, 'VariableNames', varNames);
        grpstats(democlinical_table, 'subtype', {'mean', 'std'})
        for idx = 1 : 13
            [h p ci stats] = ttest2(bayes_clinicalvar_set(group==1, idx), bayes_clinicalvar_set(group==2, idx)); 
            if(p<=0.05)
                clinicalvar_names(idx)
            end
        end
        
        FPATH = '/data3/cdb/shong/01_project/02_abcd_brain_development/03_result/';
        clinical_features = Data1;
        validation_features = Data2;
        group_ = group; % group_(group==1) = 2;  group_(group==2) = 1;  group_(group==3) = 2;  group_(group==4) = 1;
        figure;
        b = bar([ 1:5:5*16 ], ...
            [ accumarray(group_,  clinical_features(:, 1), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 2), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 3), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 4), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 5), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 6), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 7), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 8), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 9), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 10), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 11), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 12), [], @mean)'; ...
              accumarray(group_,  clinical_features(:, 13), [], @mean)'; ...
              accumarray(group_,  validation_features(:, 1), [], @mean)'; ...
              accumarray(group_,  validation_features(:, 2), [], @mean)'; ...
              accumarray(group_,  validation_features(:, 3), [], @mean)'; ], 'grouped');
        
        temp = gca; temp.XTickLabel = varNames(3:end-1); xtickangle(45);
        temp.XTick = [ 1:5:5*16 ];
        legend({'G1', 'G2', 'G3'});
        b(1).FaceColor = [0.611, 0.541, 0.878]; b(1).EdgeColor = 'none';
        b(2).FaceColor = [0.988, 0.250, 0.321]; b(2).EdgeColor = 'none';
        b(3).FaceColor = [0.988, 0.533, 0.250]; b(3).EdgeColor = 'none';
        export_fig([FPATH '/99_phenotypic_validation' ], '-m3', '-png');  close(gcf);
        
        %         'meim_exp', 'meim_com', 'meim_total',
        %         'nsc_p_ss_mean_3_items', 'fes_p_ss_fc', 'fes_y_ss_fc'
        %         'pmq_y_ss_mean', 'crpbi_y_ss_parent', 'crpbi_y_ss_caregiver'
        %         'srpf_y_ss_ses', 'srpf_y_ss_iiss', 'srpf_y_ss_dfs', 'trauma_num'
        %         'cbcl_scr_syn_internal_t', 'cbcl_scr_syn_external_t', 'cbcl_scr_syn_totprob_t'
        varNames = { 'CBCL_external_synd' };
        democlinical_var = [ ...
            bayes_validationvar_set(:, 2), ...
            group ];
        
        figure; grpstats(democlinical_var(:, 1:end-1),  group, 0.05); legend(varNames);
        
        % imaging data
        FPATH = '/data3/cdb/shong/01_project/02_abcd_brain_development/03_result/';
        imaging_feature_val =zthickness_set; range_c = [ -1 1 ];  % -1 1  |   0 3
        figure; BoSurfStatView(mean(imaging_feature_val(group==1, :), 1), ST); BoSurfStatColLim(range_c);
        figure; BoSurfStatView(mean(imaging_feature_val(group==2, :), 1), ST); BoSurfStatColLim(range_c);
        figure; BoSurfStatView(mean(imaging_feature_val(group==3, :), 1), ST); BoSurfStatColLim(range_c);
        figure; BoSurfStatView(mean(imaging_feature_val(group==4, :), 1), ST); BoSurfStatColLim(range_c);
        
        imaging_feature_val = zmyelin_set; range_c = [ -2 2 ];
        figure; BoSurfStatView(mean(imaging_feature_val(group==1, :), 1), ST); BoSurfStatColLim(range_c);
        figure; BoSurfStatView(mean(imaging_feature_val(group==2, :), 1), ST); BoSurfStatColLim(range_c);
        figure; BoSurfStatView(mean(imaging_feature_val(group==3, :), 1), ST); BoSurfStatColLim(range_c);
        figure; BoSurfStatView(mean(imaging_feature_val(group==4, :), 1), ST); BoSurfStatColLim(range_c);
        
        imaging_feature_val = bayes_myelin_set;  % bayes_thickness_set; bayes_myelin_set
        clusthres = 0.01;
        GROUP_term = term(group);
        temp = cellstr(num2str(SITE_FN));
        SITE_term = term(temp);
        M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
        M0 =  1 + AGE_term + SEX_term + SITE_term;
        slm1 = SurfStatLinMod(imaging_feature_val, M1, ST);
        slm0 = SurfStatLinMod(imaging_feature_val, M0, ST);
        slm = SurfStatF(slm1, slm0);
        [ pval, peak, clus, clusid ] = SurfStatP(slm, mask, clusthres);
        
        pmap = 1 -  fcdf(slm.t, slm.df(1), slm.df(2));
        
        figure; BoSurfStatView(slm.t, ST);
        figure; BoSurfStatView(pmap, ST, '');  colormap(mycol.red); BoSurfStatColLim([0 clusthres]);
        export_fig([FPATH '/99_thickness_ANOVA_pmap' ], '-m2', '-png');  close(gcf);
        figure; BoSurfStatView(pval.C, ST, ''); colormap(mycol.blackblue); BoSurfStatColLim([0 0.05]);
        export_fig([FPATH '/99_thickness_ANOVA_RFT' ], '-m2', '-png');  close(gcf);
        figure; BoSurfStatView(SurfStatP(slm, mask), ST);
        
        sigClus = SurfStatCluster(pval.C<0.05, ST);
        
        figure; hold on;
        mean_set = []; std_set = []; imaging_feature_set = [];
        for i = 1 : max(sigClus)
            imaging_feature_temp = mean(imaging_feature_val(:, sigClus == i), 2);
            [a b] = grpstats(imaging_feature_temp, group, { 'mean', 'std' });
            mean_set = [ mean_set a ];
            std_set     = [ std_set b ];
            imaging_feature_set = [ imaging_feature_set imaging_feature_temp ];
            eh = errorbar((1:max(group))+(i-1)*max(group), a', b');
        end
        
        feature_ind = mean(imaging_feature_set, 2);
        figure; hold on;
        for i = 1 : max(group)
            sh = scatter(i+(rand(length(feature_ind(group==i)), 1)-0.5)*0.5, feature_ind(group==i), 'filled');
            sh.SizeData = 100;
            sh.MarkerEdgeColor = 'w';
            sh.LineWidth = 1;
        end
        eh = errorbar((1:max(group)), mean(mean_set, 2)', mean(std_set, 2)');
        eh.LineStyle = 'none';
        eh.Color = [ 0 0 0 ];
        ylim([1.6 2.2]);
        ylim([0.21 0.39]);
        export_fig([RDIR '/imaging_feature_myelin_subtype_profiles' ], '-m2', '-png');  close(gcf);
        
        imaging_feature_val = cat(3, bayes_thickness_set,  bayes_myelin_set );
        clusthres = 0.01;
        
        GROUP_term = term(group);
        M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
        M0 = 1 + AGE_term + SEX_term + SITE_term;
        slm1 = SurfStatLinMod(imaging_feature_val, M1, ST);
        slm0 = SurfStatLinMod(imaging_feature_val, M0, ST);
        slm = SurfStatF(slm1, slm0);
        [ pval, peak, clus, clusid ] = SurfStatP(slm, mask);
        
        pmap = 1 -  fcdf(slm.t, slm.df(1), slm.df(2));
        
        figure; SurfStatView(slm.t, ST);
        figure; SurfStatView(pmap, ST, '');  colormap(mycol.blue); SurfStatColLim([0 clusthres]);
        figure; SurfStatView(pval.C, ST, ''); colormap(mycol.blackblue); SurfStatColLim([0 0.05]);
        figure; SurfStatView(SurfStatP(slm, mask), ST);
        
        % Phenotypic data
        pheno_feature_val = bayes_clinicalvar_set;
        slm1 = SurfStatLinMod(pheno_feature_val, M1);
        slm0 = SurfStatLinMod(pheno_feature_val, M0);
        slm = SurfStatF(slm1, slm0)
        pmap = 1 -  fcdf(slm.t, slm.df(1), slm.df(2))
        
        pheno_feature_val = bayes_validationvar_set;
        slm1 = SurfStatLinMod(pheno_feature_val, M1);
        slm0 = SurfStatLinMod(pheno_feature_val, M0);
        slm = SurfStatF(slm1, slm0)
        pmap = 1 -  fcdf(slm.t, slm.df(1), slm.df(2))
        
        pheno_feature_val = [  clinical_features(:, 3) ...
            clinical_features(:, 4) ...
            ( clinical_features(:, 5) +  clinical_features(:, 6) ) / 2 ...
            clinical_features(:, 7) ...
            ( clinical_features(:, 8) +  clinical_features(:, 9) ) / 2 ...
            ( clinical_features(:, 10) +  clinical_features(:, 11) ) / 2 ...
            clinical_features(:, 12) ...
            clinical_features(:, 13) ...
            validation_features ];
        slm1 = SurfStatLinMod(pheno_feature_val, M1);
        slm0 = SurfStatLinMod(pheno_feature_val, M0);
        slm = SurfStatF(slm1, slm0)
        pmap = 1 -  fcdf(slm.t, slm.df(1), slm.df(2))
        
    end
    
    for posthoc_prediction_analysis = 1
        
        for classification_problem = 1
           
            for SVM_multiclass_approach = 1
                
                targeted_solution = 5;
                
                if(targeted_solution == 5)
                    reordered_idx_disc = [ 1 2 3 4 5 ];
                    reordered_idx_repl = [ 1 2 3 4 5 ];
                else
                    reordered_idx_disc = [ 1 2 ];
                    reordered_idx_repl = [ 1 2 ];
                end
                
                t = templateSVM('KernelFunction', 'linear', 'OutlierFraction', 0.05);
                Mdl = fitcecoc([ imaging_features_disc1 imaging_features_disc2 ], Ci_all(1:size(imaging_features_disc1,1), targeted_solution-1), 'Learners', t);

                prediction_label = zeros(size(imaging_features_repl1, 1), 1);
                for i = 1 : size(imaging_features_repl1, 1)
                    
                    curr_case_feature1 = imaging_features_repl1(i, :);
                    curr_case_feature2 = imaging_features_repl2(i, :);
                    curr_case_feature3 = clinical_features_repl(i, :);
                    
%                     prediction_label(i) = predict(Mdl, [ curr_case_feature1 curr_case_feature2 curr_case_feature3 ]);
                    prediction_label(i) = predict(Mdl, [ curr_case_feature1 curr_case_feature2 ]);
                    
                end
                
                prediction_label_org = prediction_label;
                solution = Ci_all(size(imaging_features_disc1,1)+1:end, targeted_solution-1);
                prediction_set = sum((prediction_label-solution) == 0)/length(prediction_label)
                
                perm_iter = 100;
                prediction_rand_set = [];
                for i = 1 : perm_iter
                    i
                    randidx = randperm(size(imaging_features_disc1, 1));
                    
                    Mdl = fitcecoc([ imaging_features_disc1(randidx, :) imaging_features_disc2(randidx, :) ], Ci_all(1:size(imaging_features_disc1,1), targeted_solution-1), 'Learners', t);
                    
                    prediction_label_rand = zeros(size(imaging_features_repl1, 1), 1);
                    for i = 1 : size(imaging_features_repl1, 1)
                        
                        curr_case_feature1 = imaging_features_repl1(i, :);
                        curr_case_feature2 = imaging_features_repl2(i, :);
                        curr_case_feature3 = clinical_features_repl(i, :);
                        
                        %                     prediction_label(i) = predict(Mdl, [ curr_case_feature1 curr_case_feature2 curr_case_feature3 ]);
                        prediction_label_rand(i) = predict(Mdl, [ curr_case_feature1 curr_case_feature2 ]);
                        
                    end
                    solution = Ci_all(size(imaging_features_disc1,1)+1:end, targeted_solution-1);
                    prediction_rand_set = [ prediction_rand_set; sum((prediction_label_rand-solution) == 0)/length(prediction_label_rand) ];
                end
                
                figure; hold on;
                vs = violinplot(prediction_rand_set, ones(perm_iter, 1), 'ShowMean', true, 'ViolinColor', [0.6 0.6 0.6]); ylim([0 1]); xlim([0.5 1.5]);
                set(gcf, 'Position', [440   352   144   446]);
                s = scatter(1, prediction_set, 144, 'filled'); s.MarkerFaceColor = [1 0.5 0.5];
                export_fig(['04_subtype_classification' num2str(targeted_solution) ], '-m4', '-png');  close(gcf);
                eval([ 'prediction_rand_set_' num2str(targeted_solution) '_subtypes_solution = prediction_rand_set;']);
                
                sensitivity_class = [];
                specificity_class = [];
                for i = 1 : targeted_solution
                    
                    sensitivity_class = [ sensitivity_class sum((prediction_label == i)&(solution == i))/(sum((prediction_label == i)&(solution == i))+sum((prediction_label ~= i)&(solution == i))) ];
                    specificity_class = [ specificity_class sum((prediction_label ~= i)&(solution ~= i))/(sum((prediction_label ~= i)&(solution ~= i))+sum((prediction_label == i)&(solution ~= i))) ];
                    
                end
                
                Beta_set = [];
                for i = 1 : length(Mdl.BinaryLearners)
                   
                    Beta_set = [ Beta_set abs(Mdl.BinaryLearners{i}.Beta) ];
                    
                end
                Beta_set = mean(Beta_set, 2);
                   
                [a b] = sort(Beta_set(1:20484), 'descend');
                figure; BoSurfStatView(Beta_set(1:20484), ST);   BoSurfStatColLim([a(round(length(a)*0.2)) a(1)*0.8]); colormap([  0.8 0.8 0.8; flipud(mycol.red); ]);
                export_fig(['04_subtype_classification' num2str(targeted_solution) '_thickness_feature_map' ], '-m4', '-png');  close(gcf);
                [a b] = sort(Beta_set(20485:end), 'descend');
                figure; BoSurfStatView(Beta_set(20485:end), ST); BoSurfStatColLim([a(round(length(a)*0.2)) a(1)*0.8]); colormap([  0.8 0.8 0.8; flipud(mycol.red); ]);
                export_fig(['04_subtype_classification' num2str(targeted_solution) '_myelin_feature_map' ], '-m4', '-png');  close(gcf);
                
                temp1 = sum(Beta_set(1:20484)', 1);
                temp_yeo1 = accumarray((yeo_krienen_map+1)', temp1'); temp_yeo1(1) = [];
                temp2 = sum(Beta_set(20485:end)', 1);
                temp_yeo2 = accumarray((yeo_krienen_map+1)', temp2'); temp_yeo2(1) = [];
                
                if(curridx == 2)
                    [f, ca, o] = spider([ temp_yeo1  temp_yeo2 ],[ 'Community-wise feature contribution: externalizing' ], [0.5 2.5], ...
                        {'Visual', 'SenMot', 'DorAtt', 'Salience', 'Limbic', 'FronPar', 'DMN'}, { 'thickness', 'myelin' });
                else
                    [f, ca, o] = spider([ temp_yeo1  temp_yeo2 ],[ 'Community-wise feature contribution: externalizing' ], [0.3 1.4], ...
                        {'Visual', 'SenMot', 'DorAtt', 'Salience', 'Limbic', 'FronPar', 'DMN'}, { 'thickness', 'myelin' });
                end
                export_fig(['04_subtype_classification' num2str(targeted_solution) '_feature_map_yeo' ], '-m4', '-png');  close(gcf);
                
                
                eval( [ 'prediction_label_' num2str(targeted_solution) ' = prediction_label;' ] );
                
                thickness_pattern = zeros(curridx, 20484);
                myelin_pattern    = zeros(curridx, 20484);
                for i = 1 : targeted_solution
                    
                    thickness_pattern(i, :) = mean(imaging_features_repl1(prediction_label==i, :), 1);
                    myelin_pattern(i, :)    = mean(imaging_features_repl2(prediction_label==i, :), 1);
                    
                end
                
                if(targeted_solution == 5)
                    reordered_idx = [ 1 2 3 4 5 ];
                else
                    reordered_idx = [ 1 2 ];
                end
                
                figure; CSFSurfStatView(thickness_pattern(reordered_idx, :), ST, 'INC', [ num2str(curridx) 'Asubtypes thickness' ]); % CSFSurfStatViewColLim([1 2.6]);
                figure; CSFSurfStatView(myelin_pattern(reordered_idx, :), ST, 'INC', [ num2str(curridx) 'Asubtypes myelin' ]);  % CSFSurfStatViewColLim([0.25 0.35]);
                
                clinical_pattern   = zeros(curridx, size(clinical_features_all, 2));
                validation_pattern = zeros(curridx, size(validate_features_all, 2));
                clinicalSD_pattern   = zeros(curridx, size(clinical_features_all, 2));
                validationSD_pattern = zeros(curridx, size(validate_features_all, 2));
                
                z_clinicalvar_set_all = zscore(clinicalvar_set_all);
                z_validationvar_set_all = zscore(validationvar_set_all);
                
                legend_str = cell(1, curridx);
                for i = 1 : curridx
                    
                    legend_str{i} = [ 'subtype' num2str(i) ];
                    clinical_pattern(i, :) = mean(clinicalvar_set_all(prediction_label==i, :), 1);
                    validation_pattern(i, :) = mean(validationvar_set_all(prediction_label==i, :), 1);
                    clinicalSD_pattern(i, :) = std(clinicalvar_set_all(prediction_label==i, :), 0, 1);
                    validationSD_pattern(i, :) = std(validationvar_set_all(prediction_label==i, :), 0, 1);
                    
                end
                
                if(targeted_solution == 2)
                    
                    xoffset = 0.16;
                    
                elseif(targeted_solution == 5)
                    
                    xoffset = 0.16;
                    
                end
                
                xcoordset = round((1:targeted_solution)-(1+targeted_solution)/2)*xoffset;
                xbarcoord = [];
                for i = 1 : targeted_solution
                    
                    xbarcoord = [ xbarcoord [ 1 : size(pheno_scores_all, 2) ]' + xcoordset(i) ];
                    
                end
                
                figure; hold on;
                e = errorbar(xbarcoord, [ clinical_pattern(reordered_idx, :) validation_pattern(reordered_idx, :) ]',[ clinicalSD_pattern(reordered_idx, :) validationSD_pattern(reordered_idx, :) ]', '.');
                b = bar([ clinical_pattern(reordered_idx, :) validation_pattern(reordered_idx, :) ]'); legend(legend_str);
                for i = 1 : length(e)
                    
                    b(i).FaceColor = e(i).Color;
                    
                end
                xticks(1:length(fieldnames_set(2:end)));
                xticklabels(fieldnames_set(2:end));
                xtickangle(45);
                
            end
            
        end
        
        for subtype_info_read = 1
            
            [ subtypes_disc_k_2, text, textAndNumbers ] = xlsread('subtype_membership.xlsx', 'Discovery k=2');
            [ subtypes_repl_k_2, text, textAndNumbers ] = xlsread('subtype_membership.xlsx', 'Replication k=2');
            [ subtypes_disc_k_5, text, textAndNumbers ] = xlsread('subtype_membership.xlsx', 'Discovery k=5');
            [ subtypes_repl_k_5, text, textAndNumbers ] = xlsread('subtype_membership.xlsx', 'Replication k=5');
            
            subtypes_disc_k_2 = Ci_disc(:, 1);
            subtypes_repl_k_2 = Ci_repl(:, 1);
            subtypes_disc_k_5 = Ci_disc(:, 4);
            subtypes_repl_k_5 = Ci_repl(:, 4);
            
            thickness_disc_pattern_2 = zeros(2, 20484);
            myelin_disc_pattern_2    = zeros(2, 20484);
            thickness_repl_pattern_2 = zeros(2, 20484);
            myelin_repl_pattern_2    = zeros(2, 20484);
            thickness_disc_pattern_5 = zeros(5, 20484);
            myelin_disc_pattern_5    = zeros(5, 20484);
            thickness_repl_pattern_5 = zeros(5, 20484);
            myelin_repl_pattern_5    = zeros(5, 20484);
            for i = 1 : 2
                
                thickness_disc_pattern_2(i, :) = mean(imaging_features_disc1(subtypes_disc_k_2==i, :), 1);
                myelin_disc_pattern_2(i, :)    = mean(imaging_features_disc2(subtypes_disc_k_2==i, :), 1);
                thickness_repl_pattern_2(i, :) = mean(imaging_features_repl1(subtypes_repl_k_2==i, :), 1);
                myelin_repl_pattern_2(i, :)    = mean(imaging_features_repl2(subtypes_repl_k_2==i, :), 1);
                
            end
            
            for i = 1 : 5
                
                thickness_disc_pattern_5(i, :) = mean(imaging_features_disc1(subtypes_disc_k_5==i, :), 1);
                myelin_disc_pattern_5(i, :)    = mean(imaging_features_disc2(subtypes_disc_k_5==i, :), 1);
                thickness_repl_pattern_5(i, :) = mean(imaging_features_repl1(subtypes_repl_k_5==i, :), 1);
                myelin_repl_pattern_5(i, :)    = mean(imaging_features_repl2(subtypes_repl_k_5==i, :), 1);
                
            end
            
            reordered_idx = [ 1 2 ];
            figure; CSFSurfStatView(thickness_disc_pattern_2(reordered_idx, :), ST, 'INC', [ '2-Dsubtypes thickness' ]); % CSFSurfStatViewColLim([1 2.6]);
            figure; CSFSurfStatView(myelin_disc_pattern_2(reordered_idx, :), ST, 'INC', [ '2-Dsubtypes myelin' ]); % CSFSurfStatViewColLim([0.25 0.35]);
            
            reordered_idx = [ 1 2 ];
            figure; CSFSurfStatView(thickness_repl_pattern_2(reordered_idx, :), ST, 'INC', [ '2-Rsubtypes thickness' ]); % CSFSurfStatViewColLim([1 2.6]);
            figure; CSFSurfStatView(myelin_repl_pattern_2(reordered_idx, :), ST, 'INC', [ '2-Rsubtypes myelin' ]); % CSFSurfStatViewColLim([0.25 0.35]);
            
            reordered_idx = [ 1 2 3 4 5 ];
            figure; CSFSurfStatView(thickness_disc_pattern_5(reordered_idx, :), ST, 'INC', [ '5-Dsubtypes thickness' ]); % CSFSurfStatViewColLim([1 2.6]);
            figure; CSFSurfStatView(myelin_disc_pattern_5(reordered_idx, :), ST, 'INC', [ '5-Dsubtypes myelin' ]); % CSFSurfStatViewColLim([0.25 0.35]);
            
            reordered_idx = [ 1 2 3 4 5 ];
            figure; CSFSurfStatView(thickness_repl_pattern_5(reordered_idx, :), ST, 'INC', [ '5-Rsubtypes thickness' ]); % CSFSurfStatViewColLim([1 2.6]);
            figure; CSFSurfStatView(myelin_repl_pattern_5(reordered_idx, :), ST, 'INC', [ '5-Rsubtypes myelin' ]); % CSFSurfStatViewColLim([0.25 0.35]);
            
            for subtype2_ANCOVA = 1
                
                colormap_curr =  [ 19, 123, 194; 199, 96, 51 ]/256;
                % 2-subtype solution - discovery
                for discovery_set = 1
                    
                    SUBTYPE_discovery_str = cellstr(strcat('sub', num2str(subtypes_disc_k_2)));
                    RACE_discovery_str = cellstr(num2str(demographic_data_discovery(:, 4)));
                    SEX_discovery_str = cellstr(SEX_discovery);
                    meanthick_discovery = mean(imaging_features_disc1, 2);
                    meanmyelin_discovery = mean(imaging_features_disc2, 2);
                    
                    SEX_term = term(GENDER_FN_disc);
                    AGE_term = term(AGE_FN_disc);
                    SITE_term = term(SITE_FN_disc);
                    RACE_term = term(RACE_discovery_str);
                    GROUP_term = term(SUBTYPE_discovery_str);
                    MTHICK_term = term(meanthick_discovery);
                    MMYELIN_term = term(meanmyelin_discovery);
                    
                    thickness_yeo_ind = zeros(max(yeo_krienen_map), size(imaging_features_disc1, 1));
                    for i = 1 : max(subtypes_disc_k_2)
                       
                        thickness_yeo_ind_curr = zeros(max(yeo_krienen_map), sum(subtypes_disc_k_2==i));
                        imaging_feaures_curr = imaging_features_disc1(subtypes_disc_k_2==i, mask);
                        for j = 1 : size(imaging_feaures_curr, 1)
                            
                            temp = accumarray((yeo_krienen_map(mask)+1)', imaging_feaures_curr(j, :)', [], @mean);
                            thickness_yeo_ind_curr(:, j) = temp(2:end);
                            
                        end
                        thickness_yeo_ind(:, subtypes_disc_k_2==i) = thickness_yeo_ind_curr;
                        
                    end
                    
                    % thickness - vertex-wise
                    M1 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term;
                    slm1 = SurfStatLinMod(imaging_features_disc1, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_disc1, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2));
                    [ pval, peak, clus, clusid ] = SurfStatP( slm3, mask, 0.025 );
                    
                    resimaging_features_disc1 = imaging_features_disc1 - slm2.X(:, 2:end)*slm2.coef(2:end, :);
                    
                    [a b] = sort(subtypes_disc_k_2);
                    eta2_set = zeros(1, 20484);
                    parfor i = 1 : 20484
                        
                        i
                        
                        stats = mes1way(resimaging_features_disc1(b, i),'eta2','group',a);
                        eta2_set(i) = stats.eta2;
                        
                    end
                    
                    cohenf = sqrt(eta2_set./(1-eta2_set));
                    
                    figure; BoSurfStatView(slm3.t, ST);  colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([3 5]);
                    export_fig(['01_thickness_fstats_discovery_2_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.025]);
                    export_fig(['01_thickness_uncorrP_cohenf_discovery_2_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(pval.C, ST); BoSurfStatColLim([0 0.05]);
                    export_fig(['01_thickness_ANCOVA_pvalC_discovery_2_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(cohenf, ST); colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([0.05 0.15]);
                    export_fig(['01_thickness_ANCOVA_cohenf_discovery_2_solution' ], '-m4', '-png');  close(gcf);
                    
                    pval_thickness_subtype2_discovery = pval;
                    
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term;
                    slm2 = SurfStatLinMod(imaging_features_disc1, M2, ST);
                    res_imaging_features_disc1 = imaging_features_disc1 - slm2.X*slm2.coef;
                    
                    if(sum(pval_thickness_subtype2_discovery.C<=0.05))
                        sigver = find(pval_thickness_subtype2_discovery.C<=0.05);
                        cohenf_curr = zeros(length(sigver), 1);
                        
                        for i = 1 : length(sigver)
                            
                            currvar = res_imaging_features_disc1(:, sigver(i));
                            [a b] = sort(subtypes_disc_k_2);
                            stats = mes1way(currvar(b),'eta2','group',subtypes_disc_k_2(b));
                            eta2_set = stats.eta2;
                            cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                            
                        end
                        cohens_thickness_subtype2_discovery = zeros(1, 20484);
                        cohens_thickness_subtype2_discovery(sigver) = cohenf_curr;
                        
                    else
                        cohens_thickness_subtype2_discovery = [];
                    end
                    
                    cohenf_curr = zeros(max(yeo_krienen_map), 1);
                    yeo_feature_map = zeros(size(res_imaging_features_disc1, 1), max(yeo_krienen_map));
                    for i = 1 : max(yeo_krienen_map)
                        
                        yeo_feature_map(:, i) = sum(res_imaging_features_disc1(:, yeo_krienen_map == i), 2);
                        [a b] = sort(subtypes_disc_k_2);
                        stats = mes1way(yeo_feature_map(b, i),'eta2','group',subtypes_disc_k_2(b));
                        eta2_set = stats.eta2;
                        cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                        
                    end
                    
                    figure;
                    [f, ca, o] = spider(cohenf_curr,[ '2-subtypes solution thickness effect size (cohens f)' ], [0.0 0.05], ...
                        {'Visual', 'SenMot', 'DorAtt', 'Salience', 'Limbic', 'FronPar', 'DMN'}, {'thickness'});
                    
                    % thickness - global mean
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(meanthick_discovery, M1);
                    slm2 = SurfStatLinMod(meanthick_discovery, M2);
                    slm3  = SurfStatF(slm1, slm2); slm3.t
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2))
                    figure; vs = violinplot(meanthick_discovery, subtypes_disc_k_2, 'ShowMean', true);
                    y_lim = ylim; ylim([y_lim(1)-0.5  y_lim(2)+0.5]);
                    set(gcf, 'Position', [440   374   297   424]);
                    export_fig(['01_violin_thickness_group_mean_discovery_2_solution' ], '-m4', '-png');  close(gcf);
                    
                    [a b] = sort(subtypes_disc_k_2);
                    stats = mes1way(meanthick_discovery(b),'eta2','group',subtypes_disc_k_2(b));
                    eta2_set = stats.eta2;
                    cohenf = sqrt(eta2_set./(1-eta2_set));
                    % Cohen (1988, 285-287) proposed the following interpretation of f: f = 0.1 is a small effect, f = 0.25 is a medium effect, and f = 0.4 is a large effect.
                    
                    % thickness - rank
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(imaging_features_disc1, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_disc1, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    
                    temp1 = zeros(1, 20484);
                    temp2 = zeros(1, sum(mask));
                    [a b] = sort(slm3.t(mask), 'ascend');
                    temp2(b) = 1:sum(mask);
                    temp1(mask) = temp2;
                    figure; BoSurfStatView(temp1, ST); % colormap([ 0 0 0; flipud(mycol.red); ]); % BoSurfStatColLim([round(20484*0.5) 20484]);
                    rank_thickness_discovery_2subtypes = temp1;
                    figure; CSFSurfStatView([slm3.t; temp1], ST, 'INC', [ '2subtype-discovery' ]); % CSFSurfStatViewColLim([1 2.6]);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.025]);
                    
                    % myelin - vertex-wise
                    M1 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term;
                    slm1 = SurfStatLinMod(imaging_features_disc2, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_disc2, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2));
                    [ pval, peak, clus, clusid ] = SurfStatP( slm3, mask, 0.025 );
                    
                    resimaging_features_disc2 = imaging_features_disc2 - slm2.X(:, 2:end)*slm2.coef(2:end, :);
                    
                    [a b] = sort(subtypes_disc_k_2);
                    eta2_set = zeros(1, 20484);
                    parfor i = 1 : 20484
                        
                        i
                        
                        stats = mes1way(resimaging_features_disc2(b, i),'eta2','group',a);
                        eta2_set(i) = stats.eta2;
                        
                    end
                    
                    cohenf = sqrt(eta2_set./(1-eta2_set));
                    
                    figure; BoSurfStatView(slm3.t, ST);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.025]);
                    export_fig(['01_myelin_uncorrP_cohenf_discovery_2_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(pval.C, ST); BoSurfStatColLim([0 0.05]); colormap(mycol.blackblue);
                    export_fig(['01_myelin_ANCOVA_pvalC_discovery_2_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(cohenf, ST); colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([0.05 0.15])
                    export_fig(['01_myelin_ANCOVA_cohenf_discovery_2_solution' ], '-m4', '-png');  close(gcf);
                    
                    pval_myelin_subtype2_discovery = pval;
                    
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term;
                    slm2 = SurfStatLinMod(imaging_features_disc2, M2, ST);
                    res_imaging_features_disc2 = imaging_features_disc2 - slm2.X*slm2.coef;
                    
                    if(sum(pval_myelin_subtype2_discovery.C<=0.05))
                        
                        sigver = find(pval_myelin_subtype2_discovery.C<=0.05);
                        cohenf_curr = zeros(length(sigver), 1);
                        
                        for i = 1 : length(sigver)
                            i
                            currvar = res_imaging_features_disc2(:, sigver(i));
                            [a b] = sort(subtypes_disc_k_2);
                            stats = mes1way(currvar(b),'eta2','group',subtypes_disc_k_2(b));
                            eta2_set = stats.eta2;
                            cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                            
                        end
                        cohens_myelin_subtype2_discovery = zeros(1, 20484);
                        cohens_myelin_subtype2_discovery(sigver) = cohenf_curr;
                        
                    else
                        cohens_myelin_subtype2_discovery = [];
                    end
                    
                    cohenf_curr = zeros(max(yeo_krienen_map), 1);
                    yeo_feature_map = zeros(size(res_imaging_features_disc2, 1), max(yeo_krienen_map));
                    for i = 1 : max(yeo_krienen_map)
                        
                        yeo_feature_map(:, i) = sum(res_imaging_features_disc2(:, yeo_krienen_map == i), 2);
                        [a b] = sort(subtypes_disc_k_2);
                        stats = mes1way(yeo_feature_map(b, i),'eta2','group',subtypes_disc_k_2(b));
                        eta2_set = stats.eta2;
                        cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                        
                    end
                    
                    figure;
                    [f, ca, o] = spider(cohenf_curr,[ '2-subtypes solution myelin effect size (cohens f)' ], [0.05 0.2], ...
                        {'Visual', 'SenMot', 'DorAtt', 'Salience', 'Limbic', 'FronPar', 'DMN'}, {'myelin'});
                    
                    % myelin - global mean
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(meanmyelin_discovery, M1);
                    slm2 = SurfStatLinMod(meanmyelin_discovery, M2);
                    slm3  = SurfStatF(slm1, slm2); slm3.t
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2))
                    figure; vs = violinplot(meanmyelin_discovery, subtypes_disc_k_2, 'ShowMean', true);
                    y_lim = ylim; ylim([y_lim(1)-0.5  y_lim(2)+0.5]);
                    set(gcf, 'Position', [440   374   297   424]);
                    export_fig(['01_violin_myelin_group_mean_discovery_2_solution' ], '-m4', '-png');  close(gcf);
                    
                    % myelin - rank
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(imaging_features_disc2, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_disc2, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    
                    temp1 = zeros(1, 20484);
                    temp2 = zeros(1, sum(mask));
                    [a b] = sort(slm3.t(mask), 'ascend');
                    temp2(b) = 1:sum(mask);
                    temp1(mask) = temp2;
                    figure; BoSurfStatView(temp1, ST); % colormap([ 0 0 0; flipud(mycol.red); ]); % BoSurfStatColLim([round(20484*0.5) 20484]);
                    rank_myelin_discovery_2subtypes = temp1;
                    figure; CSFSurfStatView([slm3.t; temp1], ST, 'INC', [ '2subtype-discovery' ]); % CSFSurfStatViewColLim([1 2.6]);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.001]);
                    
                end
                
                % 2-subtype solution - replication
                for replication_set = 1
                    
                    SUBTYPE_replication_str = cellstr(strcat('sub', num2str(subtypes_repl_k_2)));
                    RACE_replication_str = cellstr(num2str(demographic_data_replication(:, 4)));
                    SEX_replication_str = cellstr(SEX_replication);
                    meanthick_replication = mean(imaging_features_repl1, 2);
                    meanmyelin_replication = mean(imaging_features_repl2, 2);
                    
                    SEX_term = term(GENDER_FN_repl);
                    AGE_term = term(AGE_FN_repl);
                    SITE_term = term(SITE_FN_repl);
                    RACE_term = term(RACE_replication_str);
                    GROUP_term = term(SUBTYPE_replication_str);
                    MTHICK_term = term(meanthick_replication);
                    MMYELIN_term = term(meanmyelin_replication);
                    
                    % thickness
                    M1 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term;
                    slm1 = SurfStatLinMod(imaging_features_repl1, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_repl1, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2));
                    [ pval, peak, clus, clusid ] = SurfStatP( slm3, mask, 0.025 );
                    
                    resimaging_features_repl1 = imaging_features_repl1 - slm2.X(:, 2:end)*slm2.coef(2:end, :);
                    
                    [a b] = sort(subtypes_repl_k_2);
                    eta2_set = zeros(1, 20484);
                    parfor i = 1 : 20484
                        
                        i
                        
                        stats = mes1way(resimaging_features_repl1(b, i),'eta2','group',a);
                        eta2_set(i) = stats.eta2;
                        
                    end
                    
                    cohenf = sqrt(eta2_set./(1-eta2_set));
                    
                    figure; BoSurfStatView(slm3.t, ST);  colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([3 5]);
                    export_fig(['01_thickness_fstats_replication_2_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.025]);
                    export_fig(['01_thickness_uncorrP_cohenf_replication_2_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(pval.C, ST); BoSurfStatColLim([0 0.05]);  colormap(mycol.blackblue);
                    export_fig(['01_thickness_ANCOVA_pvalC_replication_2_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(cohenf, ST); colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([0.05 0.2]);
                    export_fig(['01_thickness_ANCOVA_cohenf_replication_2_solution' ], '-m4', '-png');  close(gcf);
                    
                    pval_thickness_subtype2_replication = pval;
                    
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term;
                    slm2 = SurfStatLinMod(imaging_features_repl1, M2, ST);
                    res_imaging_features_repl1 = imaging_features_repl1 - slm2.X*slm2.coef;
                    
                    if(sum(pval_thickness_subtype2_replication.C<=0.05))
                        
                        sigver = find(pval_thickness_subtype2_replication.C<=0.05);
                        cohenf_curr = zeros(length(sigver), 1);
                        
                        for i = 1 : length(sigver)
                            i
                            currvar = res_imaging_features_repl1(:, sigver(i));
                            [a b] = sort(subtypes_repl_k_2);
                            stats = mes1way(currvar(b),'eta2','group',subtypes_repl_k_2(b));
                            eta2_set = stats.eta2;
                            cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                            
                        end
                        cohens_thickness_subtype2_replication = zeros(1, 20484);
                        cohens_thickness_subtype2_replication(sigver) = cohenf_curr;
                        
                    else
                        cohens_thickness_subtype2_replication = [];
                    end
                    
                    cohenf_curr = zeros(max(yeo_krienen_map), 1);
                    yeo_feature_map = zeros(size(res_imaging_features_repl1, 1), max(yeo_krienen_map));
                    for i = 1 : max(yeo_krienen_map)
                        
                        yeo_feature_map(:, i) = mean(res_imaging_features_repl1(:, yeo_krienen_map == i), 2);
                        [a b] = sort(subtypes_repl_k_2);
                        stats = mes1way(yeo_feature_map(b, i),'eta2','group',subtypes_repl_k_2(b));
                        eta2_set = stats.eta2;
                        cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                        
                    end
                    
                    figure;
                    [f, ca, o] = spider(cohenf_curr,[ '2-subtypes solution thickness effect size (cohens f)' ], [0.00 0.11], ...
                        {'Visual', 'SenMot', 'DorAtt', 'Salience', 'Limbic', 'FronPar', 'DMN'}, {'thickness'});
                    
                    % thickness - global mean
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(meanthick_replication, M1);
                    slm2 = SurfStatLinMod(meanthick_replication, M2);
                    slm3  = SurfStatF(slm1, slm2); slm3.t
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2))
                    figure; vs = violinplot(meanthick_replication, subtypes_repl_k_2, 'ShowMean', true);
                    y_lim = ylim; ylim([y_lim(1)-0.5  y_lim(2)+0.5]);
                    set(gcf, 'Position', [440   374   297   424]);
                    export_fig(['01_violin_thickness_group_mean_replication_2_solution' ], '-m4', '-png');  close(gcf);
                    
                    % thickness - rank
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(imaging_features_repl1, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_repl1, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    
                    temp1 = zeros(1, 20484);
                    temp2 = zeros(1, sum(mask));
                    [a b] = sort(slm3.t(mask), 'ascend');
                    temp2(b) = 1:sum(mask);
                    temp1(mask) = temp2;
                    figure; BoSurfStatView(temp1, ST); % colormap([ 0 0 0; flipud(mycol.red); ]); % BoSurfStatColLim([round(20484*0.5) 20484]);
                    rank_thickness_replication_2subtypes = temp1;
                    figure; CSFSurfStatView([slm3.t; temp1], ST, 'INC', [ '2subtype-replication' ]); % CSFSurfStatViewColLim([1 2.6]);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.001]);
                    
                    % myelin
                    M1 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term;
                    slm1 = SurfStatLinMod(imaging_features_repl2, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_repl2, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2));
                    [ pval, peak, clus, clusid ] = SurfStatP( slm3, mask, 0.025 );
                    
                    resimaging_features_repl2 = imaging_features_repl2 - slm2.X(:, 2:end)*slm2.coef(2:end, :);
                    
                    [a b] = sort(subtypes_repl_k_2);
                    eta2_set = zeros(1, 20484);
                    parfor i = 1 : 20484
                        
                        i
                        
                        stats = mes1way(resimaging_features_repl2(b, i),'eta2','group',a);
                        eta2_set(i) = stats.eta2;
                        
                    end
                    
                    cohenf = sqrt(eta2_set./(1-eta2_set));
                    
                    figure; BoSurfStatView(slm3.t, ST);  colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([3 5]);
                    export_fig(['01_myelin_fstats_replication_2_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.025]);
                    export_fig(['01_myelin_uncorrP_cohenf_replication_2_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(pval.C, ST); BoSurfStatColLim([0 0.05]);  colormap(mycol.blackblue);
                    export_fig(['01_myelin_ANCOVA_pvalC_replication_2_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(cohenf, ST); colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([0.05 0.2]);
                    export_fig(['01_myelin_ANCOVA_cohenf_replication_2_solution' ], '-m4', '-png');  close(gcf);
                    
                    pval_myelin_subtype2_replication = pval;
                    
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term;
                    slm2 = SurfStatLinMod(imaging_features_repl2, M2, ST);
                    res_imaging_features_repl2 = imaging_features_repl2 - slm2.X*slm2.coef;
                    
                    if(sum(pval_myelin_subtype2_replication.C<=0.05))
                        
                        sigver = find(pval_myelin_subtype2_replication.C<=0.05);
                        cohenf_curr = zeros(length(sigver), 1);
                        
                        for i = 1 : length(sigver)
                            i
                            currvar = res_imaging_features_repl2(:, sigver(i));
                            [a b] = sort(subtypes_repl_k_2);
                            stats = mes1way(currvar(b),'eta2','group',subtypes_repl_k_2(b));
                            eta2_set = stats.eta2;
                            cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                            
                        end
                        cohens_myelin_subtype2_replication = zeros(1, 20484);
                        cohens_myelin_subtype2_replication(sigver) = cohenf_curr;
                        
                    else
                        cohens_myelin_subtype2_replication = [];
                    end
                    
                    cohenf_curr = zeros(max(yeo_krienen_map), 1);
                    yeo_feature_map = zeros(size(res_imaging_features_repl2, 1), max(yeo_krienen_map));
                    for i = 1 : max(yeo_krienen_map)
                        
                        yeo_feature_map(:, i) = mean(res_imaging_features_repl2(:, yeo_krienen_map == i), 2);
                        [a b] = sort(subtypes_repl_k_2);
                        stats = mes1way(yeo_feature_map(b, i),'eta2','group',subtypes_repl_k_2(b));
                        eta2_set = stats.eta2;
                        cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                        
                    end
                    
                    figure;
                    [f, ca, o] = spider(cohenf_curr,[ '5-subtypes solution myelin effect size (cohens f)' ], [0.00 0.07], ...
                        {'Visual', 'SenMot', 'DorAtt', 'Salience', 'Limbic', 'FronPar', 'DMN'}, {'myelin'});
                    
                    % myelin - global mean
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(meanmyelin_replication, M1);
                    slm2 = SurfStatLinMod(meanmyelin_replication, M2);
                    slm3  = SurfStatF(slm1, slm2); slm3.t
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2))
                    figure; vs = violinplot(meanmyelin_replication, subtypes_repl_k_2, 'ShowMean', true);
                    y_lim = ylim; ylim([y_lim(1)-0.5  y_lim(2)+0.5]);
                    set(gcf, 'Position', [440   374   297   424]);
                    export_fig(['01_violin_myelin_group_mean_replication_2_solution' ], '-m4', '-png');  close(gcf);
                    
                    % myelin - rank
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(imaging_features_repl2, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_repl2, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    
                    temp1 = zeros(1, 20484);
                    temp2 = zeros(1, sum(mask));
                    [a b] = sort(slm3.t(mask), 'ascend');
                    temp2(b) = 1:sum(mask);
                    temp1(mask) = temp2;
                    figure; BoSurfStatView(temp1, ST); % colormap([ 0 0 0; flipud(mycol.red); ]); % BoSurfStatColLim([round(20484*0.5) 20484]);
                    rank_myelin_replication_2subtypes = temp1;
                    figure; CSFSurfStatView([slm3.t; temp1], ST, 'INC', [ '2subtype-replication' ]); % CSFSurfStatViewColLim([1 2.6]);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.001]);
                    
                end
                
                % discovery replication overlap area
                for overlap_area = 1
                    
                   temp1 = (pval_thickness_subtype2_discovery.C<=0.05 & pval_thickness_subtype2_replication.C<=0.05);
                   temp2 = (pval_myelin_subtype2_discovery.C<=0.05 & pval_myelin_subtype2_replication.C<=0.05);
                   
                   temp1_cluster = zeros(1, 20484);
                   figure; BoSurfStatView(temp1_cluster, ST); colormaptemp = colormap; colormap([0.8 0.8 0.8; colormaptemp]);
                   export_fig(['01_thickness_reproducible_mask_discovery_2_solution' ], '-m4', '-png');  close(gcf);
                   temp2_cluster = SurfStatCluster(temp2, ST);
                   figure; BoSurfStatView(temp2_cluster, ST); colormaptemp = colormap; colormap([0.8 0.8 0.8; colormaptemp]);
                   export_fig(['01_myelin_reproducible_mask_discovery_2_solution' ], '-m4', '-png');  close(gcf);
                   
                   curridx = 2;
                   for i = 1 : max(temp2_cluster)
                      
                       figure; vs = violinplot(mean(myelin_set_disc(:, temp2_cluster == i), 2), subtypes_disc_k_2, 'ShowMean', true);
                       set(gcf, 'Position', [440   374   297   424]);
                       for j = 1 : curridx
                           
                           vs(j).ViolinColor = colormap_curr(j ,:);
                           
                       end
                       
                       export_fig(['01_myelin_reproducible_mask_discovery_2_solution_profile_cluster' num2str(i) ], '-m4', '-png');  close(gcf);
                       
                   end
                   
                end
                
            end
            
            for subtype5_ANCOVA = 1
                
                colormap_curr = [ 32, 114, 176; 200, 95, 53; 237, 174, 49; 129, 77, 138; 122, 164, 68 ]/256;
                
                % 5 subtype solution - discovery
                for discovery_set = 1
                    
                    SUBTYPE_discovery_str = cellstr(strcat('sub', num2str(subtypes_disc_k_5)));
                    RACE_discovery_str = cellstr(num2str(demographic_data_discovery(:, 4)));
                    SEX_discovery_str = cellstr(SEX_discovery);
                    meanthick_discovery = mean(imaging_features_disc1, 2);
                    meanmyelin_discovery = mean(imaging_features_disc2, 2);
                    
                    SEX_term = term(GENDER_FN_disc);
                    AGE_term = term(AGE_FN_disc);
                    SITE_term = term(SITE_FN_disc);
                    RACE_term = term(RACE_discovery_str);
                    GROUP_term = term(SUBTYPE_discovery_str);
                    MTHICK_term = term(meanthick_discovery);
                    MMYELIN_term = term(meanmyelin_discovery);
                    
                    % thickness
                    M1 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term;
                    slm1 = SurfStatLinMod(imaging_features_disc1, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_disc1, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2));
                    [ pval, peak, clus, clusid ] = SurfStatP( slm3, mask, 0.025 );
                    
                    resimaging_features_disc1 = imaging_features_disc1 - slm2.X(:, 2:end)*slm2.coef(2:end, :);
                    
                    [a b] = sort(subtypes_disc_k_5);
                    eta2_set = zeros(1, 20484);
                    parfor i = 1 : 20484
                        
                        i
                        
                        stats = mes1way(resimaging_features_disc1(b, i),'eta2','group',a);
                        eta2_set(i) = stats.eta2;
                        
                    end
                    
                    cohenf = sqrt(eta2_set./(1-eta2_set));
                    
                    figure; BoSurfStatView(slm3.t, ST);  colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([3 5]);
                    export_fig(['01_thickness_fstats_discovery_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.025]);
                    export_fig(['01_thickness_uncorrP_cohenf_discovery_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(pval.C, ST); BoSurfStatColLim([0 0.05]);  colormap(mycol.blackblue);
                    export_fig(['01_thickness_ANCOVA_pvalC_discovery_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(cohenf, ST); colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([0.1 0.2]);
                    export_fig(['01_thickness_ANCOVA_cohenf_discovery_5_solution' ], '-m4', '-png');  close(gcf);
                    
                    pval_thickness_subtype5_discovery = pval;
                    
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term;
                    slm2 = SurfStatLinMod(imaging_features_disc1, M2, ST);
                    res_imaging_features_disc1 = imaging_features_disc1 - slm2.X*slm2.coef;
                    
                    if(sum(pval_thickness_subtype5_discovery.C<=0.05))
                        sigver = find(pval_thickness_subtype5_discovery.C<=0.05);
                        cohenf_curr = zeros(length(sigver), 1);
                        
                        for i = 1 : length(sigver)
                            
                            currvar = res_imaging_features_disc1(:, sigver(i));
                            [a b] = sort(subtypes_disc_k_5);
                            stats = mes1way(currvar(b),'eta2','group',subtypes_disc_k_5(b));
                            eta2_set = stats.eta2;
                            cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                            
                        end
                        cohens_thickness_subtype5_discovery = zeros(1, 20484);
                        cohens_thickness_subtype5_discovery(sigver) = cohenf_curr;
                        
                    else
                        cohens_thickness_subtype5_discovery = [];
                    end
                    
                    cohenf_curr = zeros(max(yeo_krienen_map), 1);
                    yeo_feature_map = zeros(size(res_imaging_features_disc1, 1), max(yeo_krienen_map));
                    for i = 1 : max(yeo_krienen_map)
                        
                        yeo_feature_map(:, i) = mean(res_imaging_features_disc1(:, yeo_krienen_map == i), 2);
                        [a b] = sort(subtypes_disc_k_5);
                        stats = mes1way(yeo_feature_map(b, i),'eta2','group',subtypes_disc_k_5(b));
                        eta2_set = stats.eta2;
                        cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                        
                    end
                    
                    figure;
                    [f, ca, o] = spider(cohenf_curr,[ '5-subtypes solution thickness effect size (cohens f)' ], [0.05 0.22], ...
                        {'Visual', 'SenMot', 'DorAtt', 'Salience', 'Limbic', 'FronPar', 'DMN'}, {'thickness'});
                    
                    % thickness - global mean
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(meanthick_discovery, M1);
                    slm2 = SurfStatLinMod(meanthick_discovery, M2);
                    slm3  = SurfStatF(slm1, slm2); slm3.t
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2))
                    figure; vs = violinplot(meanthick_discovery, subtypes_disc_k_5, 'ShowMean', true);
                    y_lim = ylim; ylim([y_lim(1)-0.5  y_lim(2)+0.5]);
                    set(gcf, 'Position', [440   374   297   424]);
                    export_fig(['01_violin_thickness_group_mean_discovery_5_solution' ], '-m4', '-png');  close(gcf);
                    
                    % thickness - rank
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(imaging_features_disc1, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_disc1, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    
                    temp1 = zeros(1, 20484);
                    temp2 = zeros(1, sum(mask));
                    [a b] = sort(slm3.t(mask), 'ascend');
                    temp2(b) = 1:sum(mask);
                    temp1(mask) = temp2;
                    figure; BoSurfStatView(temp1, ST); % colormap([ 0 0 0; flipud(mycol.red); ]); % BoSurfStatColLim([round(20484*0.5) 20484]);
                    rank_thickness_discovery_5subtypes = temp1;
                    figure; CSFSurfStatView([slm3.t; temp1], ST, 'INC', [ '5subtype-discovery' ]); % CSFSurfStatViewColLim([1 2.6]);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.001]);
                    
                    % myelin
                    M1 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term;
                    slm1 = SurfStatLinMod(imaging_features_disc2, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_disc2, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2));
                    [ pval, peak, clus, clusid ] = SurfStatP( slm3, mask, 0.025 );
                    
                    resimaging_features_disc2 = imaging_features_disc2 - slm2.X(:, 2:end)*slm2.coef(2:end, :);
                    
                    [a b] = sort(subtypes_disc_k_5);
                    eta2_set = zeros(1, 20484);
                    parfor i = 1 : 20484
                        
                        i
                        
                        stats = mes1way(resimaging_features_disc2(b, i),'eta2','group',a);
                        eta2_set(i) = stats.eta2;
                        
                    end
                    
                    cohenf = sqrt(eta2_set./(1-eta2_set));
                    
                    figure; BoSurfStatView(slm3.t, ST);  colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([3 5]);
                    export_fig(['01_myelin_fstats_discovery_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.025]);
                    export_fig(['01_myelin_uncorrP_cohenf_discovery_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(pval.C, ST); BoSurfStatColLim([0 0.05]);  colormap(mycol.blackblue);
                    export_fig(['01_myelin_ANCOVA_pvalC_discovery_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(cohenf, ST); colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([0.1 0.2]);
                    export_fig(['01_myelin_ANCOVA_cohenf_discovery_5_solution' ], '-m4', '-png');  close(gcf);
                    
                    pval_myelin_subtype5_discovery = pval;
                    
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term;
                    slm2 = SurfStatLinMod(imaging_features_disc2, M2, ST);
                    res_imaging_features_disc2 = imaging_features_disc2 - slm2.X*slm2.coef;
                    
                    if(sum(pval_myelin_subtype5_discovery.C<=0.05))
                        
                        sigver = find(pval_myelin_subtype5_discovery.C<=0.05);
                        cohenf_curr = zeros(length(sigver), 1);
                        
                        for i = 1 : length(sigver)
                            i
                            currvar = res_imaging_features_disc2(:, sigver(i));
                            [a b] = sort(subtypes_disc_k_5);
                            stats = mes1way(currvar(b),'eta2','group',subtypes_disc_k_5(b));
                            eta2_set = stats.eta2;
                            cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                            
                        end
                        cohens_myelin_subtype5_discovery = zeros(1, 20484);
                        cohens_myelin_subtype5_discovery(sigver) = cohenf_curr;
                        
                    else
                        cohens_myelin_subtype5_discovery = [];
                    end
                    
                    cohenf_curr = zeros(max(yeo_krienen_map), 1);
                    yeo_feature_map = zeros(size(res_imaging_features_disc2, 1), max(yeo_krienen_map));
                    for i = 1 : max(yeo_krienen_map)
                        
                        yeo_feature_map(:, i) = mean(res_imaging_features_disc2(:, yeo_krienen_map == i), 2);
                        [a b] = sort(subtypes_disc_k_5);
                        stats = mes1way(yeo_feature_map(b, i),'eta2','group',subtypes_disc_k_5(b));
                        eta2_set = stats.eta2;
                        cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                        
                    end
                    
                    figure;
                    [f, ca, o] = spider(cohenf_curr,[ '2-subtypes solution myelin effect size (cohens f)' ], [0.1 0.25], ...
                        {'Visual', 'SenMot', 'DorAtt', 'Salience', 'Limbic', 'FronPar', 'DMN'}, {'myelin'});
                    
                    % myelin - global mean
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(meanmyelin_discovery, M1);
                    slm2 = SurfStatLinMod(meanmyelin_discovery, M2);
                    slm3  = SurfStatF(slm1, slm2); slm3.t
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2))
                    figure; vs = violinplot(meanmyelin_discovery, subtypes_disc_k_5, 'ShowMean', true);
                    y_lim = ylim; ylim([y_lim(1)-0.5  y_lim(2)+0.5]);
                    set(gcf, 'Position', [440   374   297   424]);
                    export_fig(['01_violin_myelin_group_mean_discovery_5_solution' ], '-m4', '-png');  close(gcf);
                    
                    % myelin - rank
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(imaging_features_disc2, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_disc2, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    
                    temp1 = zeros(1, 20484);
                    temp2 = zeros(1, sum(mask));
                    [a b] = sort(slm3.t(mask), 'ascend');
                    temp2(b) = 1:sum(mask);
                    temp1(mask) = temp2;
                    figure; BoSurfStatView(temp1, ST); % colormap([ 0 0 0; flipud(mycol.red); ]); % BoSurfStatColLim([round(20484*0.5) 20484]);
                    rank_myelin_discovery_5subtypes = temp1;
                    figure; CSFSurfStatView([slm3.t; temp1], ST, 'INC', [ '5subtype-discovery' ]); % CSFSurfStatViewColLim([1 2.6]);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.001]);
                    
                end
                
                % 5 subtype solution - replication
                for replication_Set = 1
                    
                    SUBTYPE_replication_str = cellstr(strcat('sub', num2str(subtypes_repl_k_5)));
                    RACE_replication_str = cellstr(num2str(demographic_data_replication(:, 4)));
                    SEX_replication_str = cellstr(SEX_replication);
                    meanthick_replication = mean(imaging_features_repl1, 2);
                    meanmyelin_replication = mean(imaging_features_repl2, 2);
                    
                    SEX_term = term(GENDER_FN_repl);
                    AGE_term = term(AGE_FN_repl);
                    SITE_term = term(SITE_FN_repl);
                    RACE_term = term(RACE_replication_str);
                    GROUP_term = term(SUBTYPE_replication_str);
                    MTHICK_term = term(meanthick_replication);
                    MMYELIN_term = term(meanmyelin_replication);
                    
                    % thickness
                    M1 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term;
                    slm1 = SurfStatLinMod(imaging_features_repl1, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_repl1, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2));
                    [ pval, peak, clus, clusid ] = SurfStatP( slm3, mask, 0.025 );
                    
                    resimaging_features_repl1 = imaging_features_repl1 - slm2.X(:, 2:end)*slm2.coef(2:end, :);
                    
                    [a b] = sort(subtypes_repl_k_5);
                    eta2_set = zeros(1, 20484);
                    parfor i = 1 : 20484
                        
                        i
                        
                        stats = mes1way(resimaging_features_repl1(b, i),'eta2','group',a);
                        eta2_set(i) = stats.eta2;
                        
                    end
                    
                    cohenf = sqrt(eta2_set./(1-eta2_set));
                    
                    figure; BoSurfStatView(slm3.t, ST);  colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([3 5]);
                    export_fig(['01_thickness_fstats_replication_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.025]);
                    export_fig(['01_thickness_uncorrP_cohenf_replication_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(pval.C, ST); BoSurfStatColLim([0 0.05]);  colormap(mycol.blackblue);
                    export_fig(['01_thickness_ANCOVA_pvalC_replication_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(cohenf, ST); colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([0.1 0.2]);
                    export_fig(['01_thickness_ANCOVA_cohenf_replication_5_solution' ], '-m4', '-png');  close(gcf);
                    
                    pval_thickness_subtype5_replication = pval;
                    
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MTHICK_term;
                    slm2 = SurfStatLinMod(imaging_features_repl1, M2, ST);
                    res_imaging_features_repl1 = imaging_features_repl1 - slm2.X*slm2.coef;
                    
                    if(sum(pval_thickness_subtype5_replication.C<=0.05))
                        
                        sigver = find(pval_thickness_subtype5_replication.C<=0.05);
                        cohenf_curr = zeros(length(sigver), 1);
                        
                        for i = 1 : length(sigver)
                            i
                            currvar = res_imaging_features_repl1(:, sigver(i));
                            [a b] = sort(subtypes_repl_k_5);
                            stats = mes1way(currvar(b),'eta2','group',subtypes_repl_k_5(b));
                            eta2_set = stats.eta2;
                            cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                            
                        end
                        cohens_thickness_subtype5_replication = zeros(1, 20484);
                        cohens_thickness_subtype5_replication(sigver) = cohenf_curr;
                        
                    else
                        cohens_thickness_subtype5_replication = [];
                    end
                    
                    cohenf_curr = zeros(max(yeo_krienen_map), 1);
                    yeo_feature_map = zeros(size(res_imaging_features_repl1, 1), max(yeo_krienen_map));
                    for i = 1 : max(yeo_krienen_map)
                        
                        yeo_feature_map(:, i) = mean(res_imaging_features_repl1(:, yeo_krienen_map == i), 2);
                        [a b] = sort(subtypes_repl_k_5);
                        stats = mes1way(yeo_feature_map(b, i),'eta2','group',subtypes_repl_k_5(b));
                        eta2_set = stats.eta2;
                        cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                        
                    end
                    
                    figure;
                    [f, ca, o] = spider(cohenf_curr,[ '5-subtypes solution thickness effect size (cohens f)' ], [0.10 0.23], ...
                        {'Visual', 'SenMot', 'DorAtt', 'Salience', 'Limbic', 'FronPar', 'DMN'}, {'thickness'});
                    
                    % thickness - global mean
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(meanthick_replication, M1);
                    slm2 = SurfStatLinMod(meanthick_replication, M2);
                    slm3  = SurfStatF(slm1, slm2); slm3.t
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2))
                    figure; vs = violinplot(meanthick_replication, subtypes_repl_k_5, 'ShowMean', true);
                    y_lim = ylim; ylim([y_lim(1)-0.5  y_lim(2)+0.5]);
                    set(gcf, 'Position', [440   374   297   424]);
                    export_fig(['01_violin_thickness_group_mean_replication_5_solution' ], '-m4', '-png');  close(gcf);
                    
                    % thickness - rank
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(imaging_features_repl1, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_repl1, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    
                    temp1 = zeros(1, 20484);
                    temp2 = zeros(1, sum(mask));
                    [a b] = sort(slm3.t(mask), 'ascend');
                    temp2(b) = 1:sum(mask);
                    temp1(mask) = temp2;
                    figure; BoSurfStatView(temp1, ST); % colormap([ 0 0 0; flipud(mycol.red); ]); % BoSurfStatColLim([round(20484*0.5) 20484]);
                    rank_thickness_replication_5subtypes = temp1;
                    figure; CSFSurfStatView([slm3.t; temp1], ST, 'INC', [ '5subtype-replication' ]); % CSFSurfStatViewColLim([1 2.6]);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.001]);
                    
                    % myelin
                    M1 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term;
                    slm1 = SurfStatLinMod(imaging_features_repl2, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_repl2, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2));
                    [ pval, peak, clus, clusid ] = SurfStatP( slm3, mask, 0.025 );
                    
                    resimaging_features_repl2 = imaging_features_repl2 - slm2.X(:, 2:end)*slm2.coef(2:end, :);
                    
                    [a b] = sort(subtypes_repl_k_5);
                    eta2_set = zeros(1, 20484);
                    parfor i = 1 : 20484
                        
                        i
                        
                        stats = mes1way(resimaging_features_repl2(b, i),'eta2','group',a);
                        eta2_set(i) = stats.eta2;
                        
                    end
                    
                    cohenf = sqrt(eta2_set./(1-eta2_set));
                    
                    figure; BoSurfStatView(slm3.t, ST);  colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([3 5]);
                    export_fig(['01_myelin_fstats_replication_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.025]);
                    export_fig(['01_myelin_uncorrP_cohenf_replication_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(pval.C, ST); BoSurfStatColLim([0 0.05]);  colormap(mycol.blackblue);
                    export_fig(['01_myelin_ANCOVA_pvalC_replication_5_solution' ], '-m4', '-png');  close(gcf);
                    figure; BoSurfStatView(cohenf, ST); colormap([ 0.8 0.8 0.8; flipud(mycol.red) ]); BoSurfStatColLim([0.1 0.2]);
                    export_fig(['01_myelin_ANCOVA_cohenf_replication_5_solution' ], '-m4', '-png');  close(gcf);
                    
                    pval_myelin_subtype5_replication = pval;
                    
                    % myelin - global mean
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(meanmyelin_replication, M1);
                    slm2 = SurfStatLinMod(meanmyelin_replication, M2);
                    slm3  = SurfStatF(slm1, slm2); slm3.t
                    p = 1 - fcdf(slm3.t, slm3.df(1), slm3.df(2))
                    figure; vs = violinplot(meanmyelin_replication, subtypes_repl_k_5, 'ShowMean', true);
                    y_lim = ylim; ylim([y_lim(1)-0.5  y_lim(2)+0.5]);
                    set(gcf, 'Position', [440   374   297   424]);
                    export_fig(['01_violin_myelin_group_mean_replication_5_solution' ], '-m4', '-png');  close(gcf);
                    
                    M2 = 1 + AGE_term + SEX_term + SITE_term + MMYELIN_term;
                    slm2 = SurfStatLinMod(imaging_features_repl2, M2, ST);
                    res_imaging_features_repl2 = imaging_features_repl2 - slm2.X*slm2.coef;
                    
                    if(sum(pval_myelin_subtype5_replication.C<=0.05))
                        
                        sigver = find(pval_myelin_subtype5_replication.C<=0.05);
                        cohenf_curr = zeros(length(sigver), 1);
                        
                        for i = 1 : length(sigver)
                            i
                            currvar = res_imaging_features_repl2(:, sigver(i));
                            [a b] = sort(subtypes_repl_k_5);
                            stats = mes1way(currvar(b),'eta2','group',subtypes_repl_k_5(b));
                            eta2_set = stats.eta2;
                            cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                            
                        end
                        cohens_myelin_subtype5_replication = zeros(1, 20484);
                        cohens_myelin_subtype5_replication(sigver) = cohenf_curr;
                        
                    else
                        cohens_myelin_subtype5_replication = [];
                    end
                    
                    cohenf_curr = zeros(max(yeo_krienen_map), 1);
                    yeo_feature_map = zeros(size(res_imaging_features_repl2, 1), max(yeo_krienen_map));
                    for i = 1 : max(yeo_krienen_map)
                        
                        yeo_feature_map(:, i) = mean(res_imaging_features_repl2(:, yeo_krienen_map == i), 2);
                        [a b] = sort(subtypes_repl_k_5);
                        stats = mes1way(yeo_feature_map(b, i),'eta2','group',subtypes_repl_k_5(b));
                        eta2_set = stats.eta2;
                        cohenf_curr(i) = sqrt(eta2_set./(1-eta2_set));
                        
                    end
                    
                    figure;
                    [f, ca, o] = spider(cohenf_curr,[ '5-subtypes solution myelin effect size (cohens f)' ], [0.10 0.20], ...
                        {'Visual', 'SenMot', 'DorAtt', 'Salience', 'Limbic', 'FronPar', 'DMN'}, {'myelin'});
                    
                    % myelin - rank
                    M1 = 1 + AGE_term + SEX_term + SITE_term + GROUP_term;
                    M2 = 1 + AGE_term + SEX_term + SITE_term;
                    slm1 = SurfStatLinMod(imaging_features_repl2, M1, ST);
                    slm2 = SurfStatLinMod(imaging_features_repl2, M2, ST);
                    slm3  = SurfStatF(slm1, slm2);
                    
                    temp1 = zeros(1, 20484);
                    temp2 = zeros(1, sum(mask));
                    [a b] = sort(slm3.t(mask), 'ascend');
                    temp2(b) = 1:sum(mask);
                    temp1(mask) = temp2;
                    figure; BoSurfStatView(temp1, ST); % colormap([ 0 0 0; flipud(mycol.red); ]); % BoSurfStatColLim([round(20484*0.5) 20484]);
                    rank_myelin_replication_5subtypes = temp1;
                    figure; BoSurfStatView(p, ST); BoSurfStatColLim([0 0.001]);
                    
                end
                
                % discovery replication overlap area
                for overlap_area = 1
                    
                   temp1 = (pval_thickness_subtype5_discovery.C<=0.05 & pval_thickness_subtype5_replication.C<=0.05);
                   temp2 = (pval_myelin_subtype5_discovery.C<=0.05 & pval_myelin_subtype5_replication.C<=0.05);
                   
                   temp1_cluster = SurfStatCluster(temp1, ST);
                   figure; BoSurfStatView(temp1_cluster, ST); colormaptemp = colormap; colormap([0.8 0.8 0.8; colormaptemp]);
                   export_fig(['01_thickness_reproducible_mask_discovery_5_solution' ], '-m4', '-png');  close(gcf);
                   temp2_cluster = SurfStatCluster(temp2, ST);
                   temp2_cluster_ = temp2_cluster;
                   temp2_cluster_(temp2_cluster==1 | temp2_cluster==4 | temp2_cluster==6) = 1;
                   temp2_cluster_(temp2_cluster==7 | temp2_cluster==8) = 4;
                   temp2_cluster_(temp2_cluster==9) = 5;
                   temp2_cluster_(temp2_cluster_>5) = 0;
                   figure; BoSurfStatView(temp2_cluster_, ST); colormaptemp = colormap; colormap([0.8 0.8 0.8; colormaptemp]);
                   export_fig(['01_myelin_reproducible_mask_discovery_5_solution' ], '-m4', '-png');  close(gcf);

                  
                   curridx = 5;
                   for i = 1 : max(temp1_cluster)
                      
                       figure; vs = violinplot(mean(thickness_set_disc(:, temp1_cluster == i), 2), subtypes_disc_k_5, 'ShowMean', true);
                       set(gcf, 'Position', [440   374   297   424]);
                       for j = 1 : curridx
                           
                           vs(j).ViolinColor = colormap_curr(j ,:);
                           
                       end
                       
                       export_fig(['01_thickness_reproducible_mask_discovery_5_solution_profile_cluster' num2str(i) ], '-m4', '-png');  close(gcf);
                       
                   end
                   
                   curridx = 5;
                   for i = 1 : max(temp2_cluster_)
                      
                       figure; vs = violinplot(mean(myelin_set_disc(:, temp2_cluster == i), 2), subtypes_disc_k_5, 'ShowMean', true);
                       set(gcf, 'Position', [440   374   297   424]);
                       for j = 1 : curridx
                           
                           vs(j).ViolinColor = colormap_curr(j ,:);
                           
                       end
                       
                       export_fig(['01_myelin_reproducible_mask_discovery_5_solution_profile_cluster' num2str(i) ], '-m4', '-png');  close(gcf);
                       
                   end
                   
                end
                
            end
            
        end
        
    end
    
end