function display_and_save_matrices(graph, filename)
    % Display and save adjacency matrix
    disp('Adjacency Matrix for Graph');
    disp(graph.A);
    
    % Display and save weight matrix
    if isfield(graph, 'W')
        disp('Weight Matrix for Graph');
        disp(graph.W);
    end

    % Save matrices to a text file
    fid = fopen(filename, 'w');
    fprintf(fid, 'Adjacency Matrix for Graph:\n');
    fprintf(fid, '%g ', full(graph.A));
    fprintf(fid, '\n\nLatex:\n');
    fprintf(fid, '%s ', latex(full(graph.A) + 0 * sym('x')));
    
    if isfield(graph, 'W')
        fprintf(fid, '\n\nWeight Matrix for Graph:\n');
        fprintf(fid, '%g ', full(graph.W));
        fprintf(fid, '\n\nLatex:\n');
        fprintf(fid, '%s ', latex(full(graph.W) + 0 * sym('x')));
    end
    
    fclose(fid);
end