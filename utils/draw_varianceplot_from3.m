function draw_varianceplot_from3(para,title_string,filename,gmmhmm_projectroot,store_path,varargin)
    num_of_dists=size(para,2);
    disp(num_of_dists);
    figure('Name',['Prec. Recl. Plot for ',title_string],'Visible','On');
    hold on;
    hline = findobj(gcf, 'type', 'line');
    set(hline, 'linewidth', 3);
    styles = {'-.', '--', '-',':','-.', '--','-'};
    lengends={};
    for nn = 1:num_of_dists
        if (nargin > 5 & length(para{nn})==4)
            Distance=(1-varargin{1})*para{nn}{1}+varargin{1}*para{nn}{2};
            ground_truth_class=para{nn}{3};
            legend_string=para{nn}{4};lengends{nn}=legend_string;
        else
            Distance=para{nn}{1};
            ground_truth_class=para{nn}{2};
            legend_string=para{nn}{3};lengends{nn}=legend_string;
        end

        dist_matrix_dim=size(Distance,1);
        CLASS_NUM = 5;
        SEQ_NUM_FOR_EACH_CLASS = 10;
        dist_mean = [];
        dist_var = [];
        for jj = 1:CLASS_NUM
            third_block = (3-1)*SEQ_NUM_FOR_EACH_CLASS+1:3*SEQ_NUM_FOR_EACH_CLASS;
            jj_block = (jj-1)*SEQ_NUM_FOR_EACH_CLASS+1:jj*SEQ_NUM_FOR_EACH_CLASS;
            distances = Distance(third_block,jj_block);
            dist_mean = [ dist_mean, mean(distances(:)) ];
            dist_var = [ dist_var, std(distances(:)) ];
        end
        color=[0,0,0];
        color(mod(nn,3)+1)=0.7;
        err_bar=errorbar(dist_mean,dist_var,'color', color, 'linewidth', 3');
        hold on;
    end
	legend(lengends, 'location', 'southwest');

    %% Plot figure

    grid on;
%     title(gca,['Distance value w.r.t deviation factor of \mu'])
    title(gca,['',title_string], 'fontsize', 30);
    xlabel('deviation factor i', 'fontsize', 25);
    ylabel('Distance Value', 'fontsize', 25);
    set(gca, 'linewidth', 3, 'fontsize', 20);
    set(gca,'XTickLabel',{'-3','-2','-1','0','1','2','3'});
%     legend({'Wasserstein Distance','KL Divergence'}, 'location', 'northwest');
	legend(lengends, 'location', 'northwest');
    mkdir_if_not_exist([gmmhmm_projectroot,store_path]);
%     disp([gmmhmm_projectroot,store_path,filename,'_varianceplot.png'])
    print([gmmhmm_projectroot,store_path,filename,'_varianceplot.png'], '-dpng','-r100');
    print([gmmhmm_projectroot,store_path,filename,'_varianceplot.eps'], '-depsc','-r100');
end
