matfiles2 = ls('*.mat');

for fileindex = 1:size(matfiles2, 1)
    load(matfiles2(fileindex, :));
    trials.rtcheck = trials.rt < (mean(trials.rt) + 2* std(trials.rt));
    
    for conditionindex = 1:5
        group.performance(conditionindex, fileindex) = mean(trials.correct( trials.distance == conditionindex & trials.rtcheck ));
        group.rt(conditionindex, fileindex) = mean( trials.rt((trials.distance == conditionindex) & trials.rtcheck ));
        
    end
    
    if ~isfield(stair, 'uncrowded')
        stair.uncrowded.tilt = zeros(size(stair.crowded.tilt));
        stair.uncrowded.threshold = 0;
    end
    if ~isfield(stair.uncrowded, 'threshold')
        stair.uncrowded.threshold = mean(stair.uncrowded.tilt(size(stair.uncrowded.tilt, 2)-10:size(stair.uncrowded.tilt, 2)));
    end
    
    figure('Position', [100, 300, 1000, 700]);
    hold all
    bplt = subplot(3, 4, [1,5,9]); % tall barchart
    bar([1, 2], [stair.uncrowded.threshold, stair.crowded.threshold], 0.3);
    set(bplt, 'Ylim', [0, 45], 'XLim', [0.5, 2.5], 'XTick', [1, 2], 'XTickLabel', {'Uncrowded', 'Crowded'}, 'Box', 'off');
    title('Thresholds');
    
    perfplt = subplot(3, 4, [2,3,4,6,7,8]);
    hold all
    plot([1, 2, 3, 4, 5], [group.performance(1:5, fileindex)]);
    set(perfplt, 'Ylim', [0.4, 1], 'XLim', [0.5, 5.5], 'Box', 'off');
    title('Performance over distance');
    xlabel('Distance'); ylabel('Performance');
    
    stairplt = subplot(3, 4, [10:12]);
    hold all
    plot(abs(stair.uncrowded.tilt), '--r');
    plot(abs(stair.crowded.tilt), 'b');
    set(stairplt, 'Ylim', [0, 45], 'Box', 'off');
    
    
    group.threshold.crowded(fileindex) = stair.crowded.threshold;
    group.threshold.uncrowded(fileindex) = stair.uncrowded.threshold;
    group.staircase(fileindex, 1:size(stair.crowded.tilt, 2)) = abs(stair.crowded.tilt);
    group.personalmean(1, fileindex) = mean(trials.correct);
    group.mean.performance = mean(group.performance, 2);
    group.mean.rt = mean(group.rt, 2);
    
    suptitle(ID);
    legend(stairplt, 'Uncrowded', 'Crowded');
    
    print(gcf, '-dpng', [ID, '-crowding.png']);
end

    figure('Position', [100, 300, 1000, 700]);
    bplt = subplot(3, 4, [1,5,9]); % tall barchart
    hold all
    bar([1, 2], [mean(group.threshold.uncrowded), mean(group.threshold.crowded)], 0.3);
    scatter(ones(size(group.threshold.uncrowded)), group.threshold.uncrowded, 'ro', 'filled');
    scatter(2*ones(size(group.threshold.crowded)), group.threshold.crowded, 'ro', 'filled');
    set(bplt, 'Ylim', [0, 45], 'XLim', [0.5, 2.5], 'XTick', [1, 2], 'XTickLabel', {'Uncrowded', 'Crowded'}, 'Box', 'off');
    title('Thresholds');
    
    perfplt = subplot(3, 4, [2,3,4,6,7,8]);
    hold all
    plot([1, 2, 3, 4, 5], [group.performance(1:5, :)]);
    plot([1, 2, 3, 4, 5], [group.mean.performance], 'LineWidth', 3, 'Color', 'k');
    set(perfplt, 'Ylim', [0.5, 1], 'XLim', [0.5, 5.5], 'Box', 'off', 'XTick', [1, 2, 3, 4, 5]);
    title('Performance over distance');
    %xlabel('Distance');
    ylabel('Performance');
    
    stairplt = subplot(3, 4, [10:12]);
    hold all
    plot([1 2 3 4 5], group.rt);
    plot([1 2 3 4 5], group.mean.rt, 'LineWidth', 3, 'Color', 'k');
    set(stairplt, 'XLim', [0.5, 5.5], 'Box', 'off', 'XTick', [1 2 3 4 5]);
    ylabel('RT');
    title('Reaction Time over distance');
    
    %suptitle('Group');
    %xlabel(perfplt, 'Distance'); ylabel(perfplt, 'Performance');
    title(stairplt, 'Reaction Time over distance');
    
    title(stairplt, 'Reaction Time over distance');
    print(gcf, '-dpng', ['all-crowding.png']);
    
    close all