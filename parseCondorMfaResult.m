function [FVA_results] = parseCondorFvaResult()
    load('model.mat','model');
    load('model_thermodynamics.mat','model_thermodynamics');
    load('model_net_fluxes.mat','model_net_fluxes');    
    load('EMU.mat','EMU');
    load('idv.mat','idv');
    load('idv_known.mat','idv_known');
    load('idv_known_variance.mat','idv_known_variance');
    load('MFA.mat','MFA');

    close all;
    
    %folder name where the exection results locate
    resultDir = 'mfa_condor';
    x = CondorLoad(resultDir);
    %go over the condor running results 
    initial_random_fluxes=[];
    predicted_fluxes=[];
    predicted_concentrations=[];
    errors=[];
    for i=1:length(x)
        initial_random_fluxes       = [initial_random_fluxes x{i}.Results.initial_random_fluxes];
        predicted_fluxes            = [predicted_fluxes x{i}.Results.predicted_flux];
        predicted_concentrations    = [predicted_concentrations x{i}.Results.concentrations];
        errors                      = [errors x{i}.Results.error];
    end
    [sorted_values sorted_indexes]  = sort(errors);
    sorted_predicted_fluxes         = predicted_fluxes(:,sorted_indexes);
    sorted_predicted_concentrations = predicted_concentrations(:,sorted_indexes);
    final_predicted_flux            = sorted_predicted_fluxes(:,1);
    final_predicted_concentrattions = sorted_predicted_concentrations(:,1);
    concentration_lb = model_thermodynamics.mets_lb;
    concentration_ub = model_thermodynamics.mets_ub;
    
    
    MFA.min_model_error = sorted_values(1);
    MFA.final_predicted_flux = final_predicted_flux;
    save('MFA.mat','MFA');
    
    % plot metabolite final concentrations + lower/upper bounds
    [dR2,tR2] = xlsread('input_concentration_for_plot.xlsx');
    mets_for_plot   = tR2(1:end,1);
    mets_lb         = log10(exp(dR2(1:end,1)));
    mets_ub         = log10(exp(dR2(1:end,2)));    
    [values order]  = sort(mets_lb);
    mets_for_plot = mets_for_plot(order);
    mets_lb = mets_lb(order);
    mets_ub = mets_ub(order);
    figure;
    herrorbar((mets_lb+mets_ub)/2,[1:length(mets_for_plot)],mets_ub-(mets_lb+mets_ub)/2,'b.');
    final_predicted_concentrattions_for_plot=[];
    for(i=1:length(mets_for_plot))
        met_index = find(strcmp(mets_for_plot{i}, model_thermodynamics.mets));
        final_predicted_concentrattions_for_plot=[final_predicted_concentrattions_for_plot;log10(exp(final_predicted_concentrattions(met_index)))];
    end
    hold on;
    plot(final_predicted_concentrattions_for_plot,[1:length(mets_for_plot)],'g*','MarkerSize',15);
    hold off;
    set(gca, 'yTick', [1:length(mets_for_plot)]);
    set(gca, 'YTickLabel', mets_for_plot);
    title('Metabolite concentration bounds & predicted concentrations','FontSize',24);
    ylabel('Metabolites','FontSize',24);
    xlabel('Concentrations (log10[mM])','FontSize',24);
    set(gca, 'fontsize', 18);
    set(gca, 'XTick', [-7:1:5]);
    saveas(gcf, sprintf('./sensitivityAnalysisImg/metabolite_final_concentrations'));    
    %end of concentrations plot
    
    deltaG=[];
    RT=2.5;
    for(i=1:length(model_thermodynamics.rxns))
    %     if(~isempty(model_thermodynamics.full_rxns{i}))
        if(model_net_fluxes.is_net_flux(i))
            numerator   = sum(final_predicted_concentrattions(model_thermodynamics.product_indexes{i}));
            denominator = sum(final_predicted_concentrattions(model_thermodynamics.reactant_indexes{i}));                    
            deltaG(end+1) = model_thermodynamics.delta_G0(i)+RT*(numerator-denominator);
        else
            deltaG(end+1) = nan;
        end  
        
    end    

    title_output = sprintf('MFA');
    final_predicted_flux = MFA.final_predicted_flux;
    min_model_error = MFA.min_model_error;
    figure;
    Output(final_predicted_flux, deltaG, min_model_error, model, title_output);    
    %annotation('textarrow',x,y,'String','y = x ','LineWidth',20); %add arrow 
    saveas(gcf, sprintf('./sensitivityAnalysisImg/MFA'));
    
    
    %%
    % Plot simulated vs measured bars
    figure('units','normalized','outerposition',[0 0 1 1])
    idv_known_arr = zeros(0,1);
    for i=1:length(idv_known) 
        if ~isempty(idv_known{i})
            idv_known_arr = [idv_known_arr; i];
        end
    end
    [idv_opt idv_d cycle_error] = ComputeEmuIDV(EMU, idv, idv_known_arr, final_predicted_flux);
    % plot known metabolites - measured vs simulated
    
    for i=1:length(idv_known_arr)
        idv_opt_only=[];
        xLabelForBar=cell(0);        
        mass_isotopomer_legend=cell(0);
        x = idv_known_arr(i);
        EMU_indices = find(EMU.list(:,1)==x);
        for(j=1:max(cellfun('length',idv_known)))
            mass_isotopomer_legend{end+1}=sprintf('m+%s',num2str(j-1));
        end
%         idv_opt_only=[idv_opt_only;zeros(2,max(cellfun('length',idv_known)))];

        idv_opt_only=zeros(2,length(idv_known{x}));
        idv_opt_only(end-1,1:length(idv_known{x}))=idv_known{x};
        idv_opt_only(end,1:length(idv_opt{EMU_indices(1)}))=idv_opt{EMU_indices(1)};
        short_met_name=model.mets{x}(1:3);
        if(short_met_name=='Glu') 
            short_met_name='Gln';
        end        
        
        xLabelForBar{end+1}=strcat(short_met_name,model.mets{x}(end-2:end),'_m');
        xLabelForBar{end+1}=strcat(short_met_name,model.mets{x}(end-2:end),'_s');
               
        subplot(4,4,i);
        bar_label(idv_opt_only, idv_opt_only);
        set(gca, 'Xlim', [0.5 2.5]);
        set(gca, 'XTickLabel', xLabelForBar);
        legend(mass_isotopomer_legend, 'Location', 'NorthEastOutside');
        set(gca, 'FontSize', 11);
        set(gca, 'Ylim', [0 1.05]);
        drawnow;
    end
    suptitle('Simulated vs. Measured metabolite labeling');
    saveas(gcf, sprintf('./sensitivityAnalysisImg/simulated_vs_measured'));

    % plot unknown metabolites - simulated only
    figure('units','normalized','outerposition',[0 0 1 1]);
    index_in_subplot = 1;
    for i=1:length(model.mets)
        idv_opt_only=[];
        xLabelForBar=cell(0);    
        mass_isotopomer_legend=cell(0);
        if((ismember(i,idv_known_arr)) || (model.met_extra(i)==1))
            continue;
        end
        x = i;
        EMU_indices = find(EMU.list(:,1)==x);
        if(length(EMU_indices)==0)
            continue;
        end
        for(j=1:length(idv_opt{EMU_indices(1)}))
            mass_isotopomer_legend{end+1}=sprintf('m+%s',num2str(j-1));
        end
        idv_opt_only=zeros(2,length(idv_opt{EMU_indices(1)}));

        idv_opt_only(end-1,1:length(idv_opt{EMU_indices(1)}))=idv_opt{EMU_indices(1)};
%         idv_opt_only(end,1:length(idv_opt{EMU_indices(1)}))=idv_opt{EMU_indices(1)};
        xLabelForBar{end+1}=strcat(model.mets{x},'_s');
        xLabelForBar{end+1}=strcat(model.mets{x},'_s');
               
        subplot(3,3,index_in_subplot);
        index_in_subplot = index_in_subplot+1;
        bar_label(idv_opt_only, idv_opt_only);
        set(gca, 'Xlim', [0.5 1.5]);
        set(gca, 'XTickLabel', xLabelForBar);
        legend(mass_isotopomer_legend, 'Location', 'NorthEastOutside');
        set(gca, 'FontSize', 6);
        set(gca, 'Ylim', [0 1.05]);
        drawnow;
    end
    suptitle('Simulated metabolites');
    saveas(gcf, sprintf('./sensitivityAnalysisImg/simuilated'));
    
    
end

