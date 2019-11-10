load('model.mat','model');
load('EMU.mat','EMU');
load('idv.mat','idv');
load('idv_known.mat','idv_known');
load('idv_known_variance.mat','idv_known_variance');
load('MFA.mat','MFA');


SENSITIVITY_ANALYSIS_THRESHOLD = 3.84; %for 95% confidence interval

%folder name where the exection results locate
resultDir = 'mfa_sensitivity_analysis_condor';
x = CondorLoad(resultDir);
%go over the condor running results 
sensitivity_analysis_min_max=[];

      
for i=1:length(x)
    if((~isempty(findstr(model.rxns{i},'f'))) | (~isempty(findstr(model.rxns{i},'b'))))
%         if(~isempty(findstr(model.rxns{i},'b')))
%             continue;
%         end    
    end
    try
        min_max = [min(x{i}.Results.flux_sensitivity_analysis_vals(find(((x{i}.Results.min_model_error-MFA.min_model_error)<SENSITIVITY_ANALYSIS_THRESHOLD)==1))) max(x{i}.Results.flux_sensitivity_analysis_vals(find((((x{i}.Results.min_model_error-MFA.min_model_error)<SENSITIVITY_ANALYSIS_THRESHOLD)==1))))];
        %check lower and upper bounds and remove sensitivity analysis that
        %is outside these bounds
        if((isempty(findstr(model.rxns{i},'f'))) & (isempty(findstr(model.rxns{i},'b'))))
            min_max(1) = max(min_max(1),model.lb(i));
            min_max(2) = min(min_max(2),model.ub(i));
        end
        sensitivity_analysis_min_max = [sensitivity_analysis_min_max; min_max];        
    catch
        sensitivity_analysis_min_max = [sensitivity_analysis_min_max; NaN NaN]
    end
end 
figure;
output_sensitivity_analysis(sensitivity_analysis_min_max, MFA.min_model_error, model, 'MFA - sensitivity analysis');
saveas(gcf, sprintf('./sensitivityAnalysisImg/MFA-sensitivity analysis'));

%go over each flux and create images all fluxes, per each sensitivity
%analysis value + produce simulated vs. measured bar plot
all_possible_fluxes_array=[];
for i=1:length(x)
    try    
        x{i}.Results.flux_name;
    catch
        continue; 
    end
    flux_dir = sprintf('sensitivityAnalysisImg/%s',x{i}.Results.flux_name);
    mkdir(flux_dir);
    indexes_inside_sensitivity_analysis = find(x{i}.Results.min_model_error-MFA.min_model_error<SENSITIVITY_ANALYSIS_THRESHOLD);
    model_predicted_flux_array = x{i}.Results.model_predicted_flux_array(:,indexes_inside_sensitivity_analysis);
    all_possible_fluxes_array = [all_possible_fluxes_array model_predicted_flux_array(:,find(model_predicted_flux_array(i,:)>=model.lb(i) & model_predicted_flux_array(i,:)<=model.ub(i)))];
    min_model_error = x{i}.Results.min_model_error(indexes_inside_sensitivity_analysis);
    flux_sensitivity_analysis_vals = x{i}.Results.flux_sensitivity_analysis_vals(indexes_inside_sensitivity_analysis); 
    
    
    
    for(j=1:length(model_predicted_flux_array(1,:)))  
%     figure;
%     MFA.final_predicted_flux(j)
%     plot(model_predicted_flux_array(j,:),min_model_error);
%     hold on;
%     plot(flux_sensitivity_analysis_vals(1),min_model_error,'bo');
%     find(flux_sensitivity_analysis_vals==final_predicted_flux(1))
%     hold off;
%     %ylim([0 min_model_error*3]);
%     title(sprintf('Score/%s Flux',model.rxns{i}));
%     xlabel(sprintf('%s Flux [n-mole/(ul*h)]',model.rxns{i}));
%     ylabel('Score');
%     drawnow;
    
        
        
        
        title_output = sprintf(sprintf('MFA - %s=%.0f',x{i}.Results.flux_name,flux_sensitivity_analysis_vals(j)));
        final_predicted_flux = model_predicted_flux_array(:,j);
        final_error = min_model_error(j);

        Output(final_predicted_flux, final_error, model, title_output);  
        saveas(gcf, sprintf('./%s/%s=%.0f(MFA)',flux_dir,x{i}.Results.flux_name,flux_sensitivity_analysis_vals(j)));
        
        %%
        % Plot simulated vs measured bars
        title_output = sprintf(sprintf('simulated vs measured - %s=%.0f',x{i}.Results.flux_name,flux_sensitivity_analysis_vals(j)));
        OutputSimulatedVsMeasured(idv_known, EMU, idv, final_predicted_flux, model, title_output);
        saveas(gcf, sprintf('./%s/%s=%.0f(S&M)',flux_dir,x{i}.Results.flux_name,flux_sensitivity_analysis_vals(j)));        
        %saveas(gcf, sprintf('./sensitivityAnalysisImg/simulated_vs_measured'));

        
    end
end

save('all_possible_fluxes_array.mat','all_possible_fluxes_array');







