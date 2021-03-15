function plotPattern2D(az,el,val)
    % Plots a 2D pattern of a values as azimuth and elevation
    %
    % The parameters, `az` and `el` are the azimuth and elevation
    % angles of length naz and nel.  Value should be a matrix of size
    % `nel x naz`.
    
    % Get dimensions
    naz = length(az);
    nel = length(el);
    val = reshape(val, nel, naz);
    
    % Flip the values so that high elevations are on the top
    val = flipud(val);
    el = flipud(el);
    
    % Plot the result
    imagesc(az, el, val);
    colorbar();
    
    
    