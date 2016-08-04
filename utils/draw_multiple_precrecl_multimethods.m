function draw_multiple_precrecl_multimethods(para,title_string,filename,gmmhmm_projectroot,store_path,varargin)
    num_of_class=size(para,2);
    num_of_precrecl_lines=size(para{1},2);
    disp(num_of_precrecl_lines);
    figure('Name',['Prec. Recl. Plot for',title_string],'Visible','On');
    hold on;
    hline = findobj(gcf, 'type', 'line');
    set(hline, 'linewidth', 3);
    styles = {'-.', '--', '-',':','-.', '--','-'};
    lengends={};
    for class_idx = 1:num_of_class
        for nn = 1:num_of_precrecl_lines
            if (nargin > 5 & length(para{class_idx}{nn})==4)
                disp('hello');
                Distance=(1-varargin{1})*para{class_idx}{nn}{1}+varargin{1}*para{class_idx}{nn}{2};
                ground_truth_class=para{class_idx}{nn}{3};
                legend_string=para{class_idx}{nn}{4};lengends{(class_idx-1)*num_of_precrecl_lines+nn}=legend_string;
            else
                Distance=para{class_idx}{nn}{1};
                ground_truth_class=para{class_idx}{nn}{2};
                legend_string=para{class_idx}{nn}{3};lengends{(class_idx-1)*num_of_precrecl_lines+nn}=legend_string;
            end

            dist_matrix_dim=size(Distance,1);

            prec=zeros(dist_matrix_dim,size(1:max( round(dist_matrix_dim / 100), 1):dist_matrix_dim,2));
            recl=zeros(dist_matrix_dim,size(1:max( round(dist_matrix_dim / 100), 1):dist_matrix_dim,2));

            for i=1:dist_matrix_dim
                Distance(i,i)=0;
                dist_to_others=-Distance(i,:);
                dist_to_others(i)=-Inf;
                ground_truth_label_i=zeros(1,dist_matrix_dim);
                ground_truth_label_i(ground_truth_class==ground_truth_class(i))=1;
                [~,~, prec(i,:), recl(i,:), ~, ~] = precisionRecall(dist_to_others,ground_truth_label_i);
            end
            color=[0,0,0];
            color(mod(nn,3)+1)=0.7;
%             [0,20*nn,30+50*nn]/255
            plot(mean(recl(:,:),1),mean(prec(:,:),1), 'color', color, 'linewidth', 3','linestyle',styles{class_idx})  

        end
    end
    box on;
    grid on;
    
    title(gca,['Prec. Recl. Plot for ',title_string]);
    xlabel('Recall', 'fontsize', 20);
    ylim([0,1]);
    ylabel('Precision', 'fontsize', 20);
    set(gca, 'linewidth', 3, 'fontsize', 15);
    disp(lengends);
	legend(lengends, 'location', 'northeast');%'southwest'
    mkdir_if_not_exist([gmmhmm_projectroot,store_path]);
    print([gmmhmm_projectroot,store_path,filename,'_prec_recl_overall.png'], '-dpng','-r300');
    print([gmmhmm_projectroot,store_path,filename,'_prec_recl_overall.eps'], '-depsc','-r300');
end
