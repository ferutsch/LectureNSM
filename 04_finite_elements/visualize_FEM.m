function visualize_temp(nodalCoords, edof, result)

    figure;
    hold on;
    %colormap(jet);
    colorVals = result(:); % Ensure column vector

    % Loop through each element and plot
    for el = 1:size(edof, 1)
        nodeIndices = edof(el, :);
        coords = nodalCoords(nodeIndices, :);
        vals = colorVals(nodeIndices);

        % Close the polygon by repeating the first node
        coords = [coords; coords(1,:)];
        vals = [vals; vals(1)];

        % Patch for colored element
        patch('XData', coords(:,1), ...
              'YData', coords(:,2), ...
              'FaceVertexCData', vals, ...
              'FaceColor', 'interp', ...
              'EdgeColor', [0.3 0.3 0.3]); % optional: gray mesh lines
    end

    % Plot nodal points
    plot(nodalCoords(:,1), nodalCoords(:,2), 'k.', 'MarkerSize', 8);

    % Enhance plot
    axis equal
    xlim ([-0.5 0.5])
    ylim ([-0.5 0.5])
    colorbar;
    title('Temperature Distribution');
    xlabel('X');
    ylabel('Y');
    grid on;
end