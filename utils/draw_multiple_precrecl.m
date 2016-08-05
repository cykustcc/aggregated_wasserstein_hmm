function draw_multiple_precrecl(para,title_string,filename,gmmhmm_projectroot,store_path,varargin)
    num_of_precrecl_lines=size(para,2);
%     disp(num_of_precrecl_lines);
    figure('Name',['Prec. Recl. Plot for ',title_string],'Visible','On');
    hold on;
    hline = findobj(gcf, 'type', 'line');
    set(hline, 'linewidth', 3);
    styles = {'-.', '--', '-',':','-.', '--','-'};
    lengends={};
    for nn = 1:num_of_precrecl_lines
        if (nargin > 5 & length(para{nn})==4)
%             disp('hello');
            Distance=(1-varargin{1})*para{nn}{1}+varargin{1}*para{nn}{2};
            ground_truth_class=para{nn}{3};
            legend_string=para{nn}{4};lengends{nn}=legend_string;
        else
            Distance=para{nn}{1};
            ground_truth_class=para{nn}{2};
            legend_string=para{nn}{3};lengends{nn}=legend_string;
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
%             disp('dist_to_others:');
%             disp(length(dist_to_others));
%             disp('ground_truth_label_i:');
%             disp(length(ground_truth_label_i));
            [~,~, prec(i,:), recl(i,:), ~, ~] = precisionRecall(dist_to_others,ground_truth_label_i);
        end
        color=[0,0,0];
        color(mod(nn,3)+1)=0.7;
        plot(mean(recl(:,:),1),mean(prec(:,:),1), 'color', color, 'linewidth', 3','linestyle',styles{mod(nn,3)+1})  
    
    end
    box on;
    grid on;
    
%     title(gca,['Prec. Recl. Plot for ',title_string]);
    title(gca,['',title_string], 'fontsize', 30);
    xlabel('Recall', 'fontsize', 20);
    ylim([0,1]);
    ylabel('Precision', 'fontsize', 20);
    set(gca, 'linewidth', 3, 'fontsize', 20);
%     disp(lengends);
	legend(lengends, 'location', 'southwest');
    mkdir_if_not_exist([gmmhmm_projectroot,store_path]);
    print([gmmhmm_projectroot,store_path,filename,'_prec_recl_overall.png'], '-dpng','-r300');
    print([gmmhmm_projectroot,store_path,filename,'_prec_recl_overall.eps'], '-depsc','-r300');
end



