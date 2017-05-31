function [patients_sigma_pred,controls_sigma_pred,patients_threshold,controls_threshold] = optimiseperception_threshold_forpublication
%Must run 'Readdemographics.m' first to get the right inputs

global headings patientdata controldata mu_pred sigma_pred mu_input sigma_input threshold
global vocode_report_4_patient vocode_report_8_patient vocode_report_16_patient vocode_report_4_control vocode_report_8_control vocode_report_16_control normalised_model_meansarray 

% Read in the clarity of input based on vocode report task
vocode_report_4_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'Vocode Report 4'))));
vocode_report_8_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'Vocode Report 8'))));
vocode_report_16_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'Vocode Report 16'))));


vocode_report_4_control = controldata(:,(~cellfun('isempty',strfind(headings,'Vocode Report 4'))));
vocode_report_8_control = controldata(:,(~cellfun('isempty',strfind(headings,'Vocode Report 8'))));
vocode_report_16_control = controldata(:,(~cellfun('isempty',strfind(headings,'Vocode Report 16'))));

% Read in the clarity ratings from the MEG task (and neutral task, in case
% I analyse this later)

match_rating_4_chan_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'MEG Match4'))));
match_rating_8_chan_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'MEG Match8'))));
match_rating_16_chan_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'MEG Match16'))));

mismatch_rating_4_chan_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'MEG Mismatch4'))));
mismatch_rating_8_chan_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'MEG Mismatch8'))));
mismatch_rating_16_chan_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'MEG Mismatch16'))));

neutral_rating_4_chan_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'MEG Neutral4'))));
neutral_rating_8_chan_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'MEG Neutral8'))));
neutral_rating_16_chan_patient = patientdata(:,(~cellfun('isempty',strfind(headings,'MEG Neutral16'))));


match_rating_4_chan_control = controldata(:,(~cellfun('isempty',strfind(headings,'MEG Match4'))));
match_rating_8_chan_control = controldata(:,(~cellfun('isempty',strfind(headings,'MEG Match8'))));
match_rating_16_chan_control = controldata(:,(~cellfun('isempty',strfind(headings,'MEG Match16'))));

mismatch_rating_4_chan_control = controldata(:,(~cellfun('isempty',strfind(headings,'MEG Mismatch4'))));
mismatch_rating_8_chan_control = controldata(:,(~cellfun('isempty',strfind(headings,'MEG Mismatch8'))));
mismatch_rating_16_chan_control = controldata(:,(~cellfun('isempty',strfind(headings,'MEG Mismatch16'))));

neutral_rating_4_chan_control = controldata(:,(~cellfun('isempty',strfind(headings,'MEG Neutral4'))));
neutral_rating_8_chan_control = controldata(:,(~cellfun('isempty',strfind(headings,'MEG Neutral8'))));
neutral_rating_16_chan_control = controldata(:,(~cellfun('isempty',strfind(headings,'MEG Neutral16'))));

% Work out the priors
mu_pred = 10; %Arbitrarily place prediction and input on x-scale. When 'match', mu_input will be the same, when mismatch, mu_input will be different

global connum patnum mu_input_match mu_input_mismatch mu_pred meansarray
connum = 0; %for patternsearch

for patnum = 1:size(patientdata,1)
    %     for patnum = 3
    % First find optimal fit arguments
    meansarray = [mismatch_rating_4_chan_patient(patnum), mismatch_rating_8_chan_patient(patnum), mismatch_rating_16_chan_patient(patnum), neutral_rating_4_chan_patient(patnum), neutral_rating_8_chan_patient(patnum), neutral_rating_16_chan_patient(patnum), match_rating_4_chan_patient(patnum), match_rating_8_chan_patient(patnum), match_rating_16_chan_patient(patnum)];
    denormed_meansarray = (meansarray-1)./3; %De-normalise
    mu_input_match = 10;
    mu_input_mismatch = 20;
 
    [model_arguments, fval] = patternsearch(@search_subfunction_bayes_modelling_sohoglu2016_forpublication,[3,0],[],[],[],[],[eps, 0],[10, 1]);

    sigma_pred = model_arguments(1);
    threshold = model_arguments(2);
    
    %Now model MEG Clarity ratings (normalised_model_meansarray is a
    %global, from the optimisation function)
    eval(['modelled_error_for_patient' num2str(patnum) ' = search_subfunction_bayes_modelling_sohoglu2016_forpublication(model_arguments)'])
    all_meansarray = [normalised_model_meansarray,meansarray];
    stesarray = zeros(size(all_meansarray)); %Dummy stes for now.
        figure
        set(gcf,'position',[100,100,1200,800])
        barweb(all_meansarray,stesarray,[],{'4 channels';'8 channels';'16 channels'},['Compared Clarity Ratings by Prime Type and Vocoder Channels for Patient ' num2str(patnum)],[],'Mean Clarity Rating',[],[],{'Model_Match','Model_Mismatch','Actual_Match','Actual_Mismatch'}) ;
        legend('Model Match','Model Mismatch','Actual Match','Actual Mismatch','location','NorthWest');
        set(gca,'ylim',[1,4]);
    
    all_meansarray_patients(patnum,:,:) = all_meansarray;
    
    patients_sigma_pred(patnum) = sigma_pred;
    patients_threshold(patnum) = threshold;
end
patnum = 0; %for patternsearch
for connum = 1:size(controldata,1)
    meansarray = [mismatch_rating_4_chan_control(connum), mismatch_rating_8_chan_control(connum), mismatch_rating_16_chan_control(connum), neutral_rating_4_chan_control(connum), neutral_rating_8_chan_control(connum), neutral_rating_16_chan_control(connum), match_rating_4_chan_control(connum), match_rating_8_chan_control(connum), match_rating_16_chan_control(connum)];
    denormed_meansarray = (meansarray-1)./3; %De-normalise
    mu_input_match = 10;
    mu_input_mismatch = 20;
    [model_arguments, fval] = patternsearch(@search_subfunction_bayes_modelling_sohoglu2016_forpublication,[3,0],[],[],[],[],[eps, 0],[10, 1]);
    sigma_pred = model_arguments(1);
    threshold = model_arguments(2);
    
    %Now model MEG Clarity ratings
    eval(['modelled_error_for_control' num2str(connum) ' = search_subfunction_bayes_modelling_sohoglu2016_forpublication(model_arguments)'])
    
    all_meansarray = [normalised_model_meansarray,meansarray];
    
    stesarray = zeros(size(all_meansarray)); %Dummy stes for now.
        figure
        set(gcf,'position',[100,100,1200,800])
        barweb(all_meansarray,stesarray,[],{'4 channels';'8 channels';'16 channels'},['Compared Clarity Ratings by Prime Type and Vocoder Channels for control ' num2str(connum)],[],'Mean Clarity Rating',[],[],{'Model_Match','Model_Mismatch','Actual_Match','Actual_Mismatch'}) ;
        legend('Model Match','Model Mismatch','Actual Match','Actual Mismatch','location','NorthWest');
        set(gca,'ylim',[1,4]);
    
    all_meansarray_controls(connum,:,:) = all_meansarray;
    
    controls_sigma_pred(connum) = sigma_pred;
    controls_threshold(connum) = threshold;
end
figure
stesarray = reshape(stesarray,3,6);
barweb(reshape(mean(all_meansarray_patients,1),3,6),stesarray,[],{'4 channels';'8 channels';'16 channels'},['Compared Clarity Ratings by Prime Type and Vocoder Channels for All Patients'],[],'Mean Clarity Rating',[],[],{}) ;
set(gca,'ylim',[1,4]);
figure
barweb(reshape(mean(all_meansarray_controls,1),3,6),stesarray,[],{'4 channels';'8 channels';'16 channels'},['Compared Clarity Ratings by Prime Type and Vocoder Channels for All Controls'],[],'Mean Clarity Rating',[],[],{}) ;
set(gca,'ylim',[1,4]);

patient_means_forline = reshape(mean(all_meansarray_patients,1),3,6);
control_means_forline = reshape(mean(all_meansarray_controls,1),3,6);

figure
set(gcf,'position',[100,100,1200,800])
lineplot = tight_subplot(1,2,[0 0],[.1 .1],[.1 .1]);
axes(lineplot(1));
set(gca,'ylim',[1,4])
errorbar(control_means_forline(:,1),stesarray(:,1),'r--','linewidth',3);
hold on
errorbar(control_means_forline(:,3),stesarray(:,3),'k--','linewidth',3);
errorbar(control_means_forline(:,4),stesarray(:,4),'r','linewidth',3);
%errorbar(control_means_forline(:,5),stesarray(:,5),'g','linewidth',3);
errorbar(control_means_forline(:,6),stesarray(:,6),'k','linewidth',3);
set(gca,'ylim',[1,4],'LineWidth', 2, 'Xtick', [1 2 3], 'XTickLabel',[4,8,16],'Fontsize',[14],'FontName','Tahoma','YAxisLocation','right')
secondlegend = legend('Model MisMatch','Model Match','MisMatch','Match','location','SouthEast');
set(secondlegend,'FontSize',18);
title('Controls','Color','k','fontsize',20)
ylabel('Clarity Rating')
xlabel('Vocode Channels')

axes(lineplot(2));
set(gca,'ylim',[1,4])
errorbar(patient_means_forline(:,1),stesarray(:,1),'r--','linewidth',3);
hold on
errorbar(patient_means_forline(:,3),stesarray(:,3),'k--','linewidth',3);
errorbar(patient_means_forline(:,4),stesarray(:,4),'r','linewidth',3);
errorbar(patient_means_forline(:,6),stesarray(:,6),'k','linewidth',3);
set(gca,'ylim',[1,4],'LineWidth', 2, 'Xtick', [1 2 3], 'XTickLabel',[4,8,16],'Fontsize',[14],'FontName','Tahoma','YAxisLocation','right')
set(secondlegend,'FontSize',18);
title('Patients','Color','k','fontsize',20)
ylabel('Clarity Rating')
xlabel('Vocode Channels')
img = getframe(gcf);
try %This is an absolute path, so clearly won't work anywhere except my laptop
addpath('C:\Users\Thoma\Documents\Academic work\MATLAB\ojwoodford-export_fig-216b30e')
export_fig Model_Fits_NoNeutral.png -transparent
export_fig Model_Fits_NoNeutral.pdf -transparent
catch
end


figure
boxplot([patients_sigma_pred',controls_sigma_pred'],'plotstyle','compact','symbol','k.','medianstyle','line','colors',hsv2rgb([0 0.6 0.6; 0.3 0.6 0.6; 0.6 0.6 0.6]),'jitter',0,'widths',0.5)
title('Standard Deviation of Prior')
xtix = {'\bfnfvPPA','\bfControls'};   % Your labels
xtixloc = [1 2];      % Your label locations
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc,'FontSize',15,'ylim',[0 3]);
set(findobj(gca,'Type','text'),'FontSize',80);
set(findobj(gcf,'Tag','Median'),'Color',[0 0 0],'LineWidth',4);
set(findobj(gcf,'Tag','Outliers'),'MarkerEdgeColor',[0 0 0],'MarkerSize',30);
set(findobj(gcf,'Tag','Box'),'linewidth',20);
set(findobj(gcf,'Tag','Whisker'),'linewidth',4);
ylabel('\bfA.U.')
try %This is an absolute path, so clearly won't work anywhere except my laptop
    tc_sigstar([1,2],ranksum(patients_sigma_pred,controls_sigma_pred))
    export_fig Standard_Deviation_of_Prior.png -transparent
    export_fig Standard_Deviation_of_Prior.pdf -transparent
catch
end

figure
boxplot([patients_threshold',controls_threshold'],'plotstyle','compact','symbol','k.','medianstyle','line','colors',hsv2rgb([0 0.6 0.6; 0.3 0.6 0.6; 0.6 0.6 0.6]),'jitter',0,'widths',0.5)
title('Perceptual Threshold')
xtix = {'\bfnfvPPA','\bfControls'};   % Your labels
xtixloc = [1 2];      % Your label locations
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc,'FontSize',15,'ylim',[0 0.4]);
set(findobj(gca,'Type','text'),'FontSize',80);
set(findobj(gcf,'Tag','Median'),'Color',[0 0 0],'LineWidth',4);
set(findobj(gcf,'Tag','Outliers'),'MarkerEdgeColor',[0 0 0],'MarkerSize',30);
set(findobj(gcf,'Tag','Box'),'linewidth',20);
set(findobj(gcf,'Tag','Whisker'),'linewidth',4);
ylabel('\bfA.U.')
try %This is an absolute path, so clearly won't work anywhere except my laptop
export_fig Perceptual_Threshold.png -transparent
export_fig Perceptual_Threshold.pdf -transparent
catch
end

nonpara_pval_sigma_pred = ranksum(patients_sigma_pred,controls_sigma_pred)
nonpara_pval_threshold = ranksum(patients_threshold,controls_threshold)
nonpara_pval_ratios = ranksum(patients_sigma_pred./patients_threshold,controls_sigma_pred./controls_threshold)

[~, para_pval_sigma_pred] = ttest2(patients_sigma_pred,controls_sigma_pred,'vartype','unequal')
[~, para_pval_threshold] = ttest2(patients_threshold,controls_threshold,'vartype','unequal')

