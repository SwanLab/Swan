function plotRmseOverEpisodes(rmse_vector)
    figure;
    plot(1:length(rmse_vector), rmse_vector, 'b-', 'LineWidth', 1.5);
    xlabel('Episode');
    ylabel('RMSE to Reference Q');
    title('RMSE of Estimated Q Compared to Reference');
    grid on;
end