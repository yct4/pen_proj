clear all;
close all;
data_sheet_path = 'C:\Users\nico\pen_proj\';
folder = '';

%folder = 'test_0705\';
file = 'letter_test_0705'; % abcde
file = 'letter_test_0705_2'; % abcde
file = 'test1';


path = [data_sheet_path folder file '.csv'];
process_file(path);

function process_file(path)
    q_in = [];
    g_input = readtable(path);

    M = table2array(g_input);
    idx = all(isnan(M),2);
    idr = diff(find([1;diff(idx);1]));
    D = mat2cell(M,idr(:),size(M,2));
    D{1:2:end};

    j = 0;

    for m = 1:size(D)
        input_arr = D(m);
        in = cell2mat(input_arr);

        if size(in, 1) < 20 % TODO threshold, was 10
            continue;
        end
        j = j + 1;
%         if (j == 1) % skip first letter bc not part of recording
%             continue;
%         end
        
        quat_in = in(1:end,1:end) / 1073741824;
        q1 = quat_in(1:end,1);
        q2 = quat_in(1:end,2);
        q3 = quat_in(1:end,3);
        q0 = sqrt(abs(1.0 - ((q1.^2) + (q2.^2) + (q3.^2))));
        quat_in = [q0 q1 q2 q3];
        
%         print = ['num = ' num2str(j) '; len = ' num2str(size(q1,1))] % letter number and length
        
        process_data(quat_in, j);
        q_in = [q_in; q0 q1 q2 q3];
%         break; % one letter only
    end
%     process_data(q_in, j); % whole word
end

function process_data(q, num)
    srt_i = 1;
    end_i = size(q(:,1),1);

    qw = q(:,1);
    qx = q(:,2);
    qy = q(:,3);
    qz = q(:,4);

    angle = 2 * acos(qw);
    x = qx ./ sqrt(1-qw.^2);
    y = qy ./ sqrt(1-qw.^2);
    z = qz ./ sqrt(1-qw.^2);
    
    x_i = x .* angle;
    y_i = y .* angle;
    z = z .* angle;

    y = x_i;
    x = y_i;
    
    % INITIAL SCALING
    a = 1;
    b = 5;
    [x, y] = scale_xy(x, y, a, b);
    z_scaled = (b - a) .* (z - min(z)) ./ (max(z) - min(z)) + a;
    z = z_scaled;
    
    % SMOOTHING
    degree = 5; % THRESHOLD moving avg
    x_fit = smooth(x, degree);
    y_fit = smooth(y, degree);
    z_fit = smooth(z, degree);
    
    % PLOTTING
%     figure; hold on; title(['x,y raw and smooth. letter: ' num2str(num)]); grid on;
%     plot(x_fit, '-b', 'LineWidth', 2);
%     plot(y_fit, '-r', 'LineWidth', 2);
%     plot(z_fit, '-g', 'LineWidth', 2);
% %     yyaxis right; grid on;
%     legend('x','y','z');
%     hold off;
    
%     figure; hold on; grid on; title(['xyz b4 rotation. letter: ' num2str(num)]);
%     view([0 20]);
%     xlabel('x-axis');
%     ylabel('y-axis');
%     zlabel('z-axis');
%     plot3(x_fit(srt_i:end_i),y_fit(srt_i:end_i),z_fit(srt_i:end_i));
%     hold off;

    % ROTATION y-AXIS
    R = roty(-45);
    Ar_mat = [];
    for i = 1:1:size(x_fit,1)
        Ab(1,1) = x_fit(i);
        Ab(2,1) = y_fit(i);
        Ab(3,1) = z_fit(i);
        Ar = (R*Ab).';
        Ar_mat = [Ar_mat; Ar];  
    end
    
    x_fit = Ar_mat(:,1);
    y_fit = Ar_mat(:,2);
    z_fit = Ar_mat(:,3);

    x_fit = -x_fit;
    [x_fit, y_fit] = scale_xy(x_fit, y_fit, a, b);
    z_scaled = (b - a) .* (z_fit - min(z_fit)) ./ (max(z_fit) - min(z_fit)) + a;
    z_fit = z_scaled;
    
    figure; hold on; title(['x,y raw and smooth. letter: ' num2str(num)]); grid on;
%     plot(z_fit(srt_i:end_i),x_fit(srt_i:end_i), '-b', 'LineWidth', 2); % original, no truncation
    
    % truncate at first z-maximum
    zmax = islocalmax(z_fit);
    zmax_indices = find(zmax);
    
    if (size(zmax_indices,1) > 0)
        start_sample = zmax_indices(1);
    else 
        start_sample = 1;
    end
    
    x_fit = x_fit(start_sample:end); % truncate
    y_fit = y_fit(start_sample:end); % truncate
    z_fit = z_fit(start_sample:end); % truncate
    [x_fit, y_fit] = scale_xy(x_fit, y_fit, a, b);
    z_scaled = (b - a) .* (z_fit - min(z_fit)) ./ (max(z_fit) - min(z_fit)) + a;
    z_fit = z_scaled;
    
    plot(z_fit,x_fit, '-r', 'LineWidth', 2);
    
    % if first z value is 5, truncate 
    zmax = islocalmax(z_fit);
    zmax_indices = find(zmax);
    zmax_vals = z_fit(zmax_indices);
    cmp = zmax_vals < 5;
    temp = find(cmp);
    th = max(zmax_vals(temp));
    
    if (z_fit(1) == b && size(th,1) > 0)
        cmp = z_fit < th;
        temp = find(cmp);
        start_sample = temp(1);
        x_fit = x_fit(start_sample:end); % truncate
        y_fit = y_fit(start_sample:end); % truncate
        z_fit = z_fit(start_sample:end); % truncate
        [x_fit, y_fit] = scale_xy(x_fit, y_fit, a, b);
        z_scaled = (b - a) .* (z_fit - min(z_fit)) ./ (max(z_fit) - min(z_fit)) + a;
        z_fit = z_scaled;
    end
    
    plot(z_fit,x_fit, '-og', 'LineWidth', 2);
    
    [z_fit,x_fit] = normalize_char(z_fit,x_fit,num);
    zx_dist = calc_dist(z_fit,x_fit)

    p1 = size(x_fit,1);
    resample1 = ['num = ' num2str(num), '; len = ' num2str(p1)]
    plot(z_fit,x_fit,'o-', 'LineWidth', 2); 
    hold off;
%     zx_dist = calc_dist(z_fit,x_fit) % distance isnt perfectly uniformly spaced but mostly uniform

    chaincode = chain_code(z_fit,x_fit)
    

end

function chaincode = chain_code(x,y,num) % parameterizable by # of directions
    num_dir = 8; % number of directions, assumes divisible by 4 
    slice_angle = 2*pi / num_dir;

    x_dist = x(2:end) - x(1:end-1);
    y_dist = y(2:end) - y(1:end-1);
    angle = atan(x_dist./y_dist)
    chaincode = zeros(size(angle,1),1)-1;
    
    % upper half
    cmp = y_dist >= 0;
    temp = find(cmp);
    chaincode(temp) = mod(floor(angle(temp) / slice_angle),8);
    
    % lower half
    cmp = y_dist < 0;
    temp = find(cmp);
    chaincode(temp) = mod(4 + floor(angle(temp) / slice_angle),8);
    
%     size(chaincode)
end

function [ret_x,ret_y] = normalize_char(x,y,num)
    resample_len = 30; % THRESHOLD
    
    xy_dist = calc_dist(x,y);
    pt_dist = sum(xy_dist) / (resample_len);

    [z_norm,x_norm] = resample_char(x, y, pt_dist, resample_len);
    plot(z_norm,x_norm, '-om', 'LineWidth', 1);
    
    z_temp = spline(1:1:size(z_norm,1), z_norm);
    x_temp = spline(1:1:size(x_norm,1), x_norm);
    xx = 0:size(z_norm,1)/(resample_len):size(z_norm,1);

    z_rs = ppval(z_temp,xx);
    z_rs = z_rs(2:end).';
    x_rs = ppval(x_temp,xx);
    x_rs = x_rs(2:end).';
    
    ret_x = z_rs;
    ret_y = x_rs;
    
end

function xy_dist = calc_dist(x,y)
    x_dist = abs(x(2:end) - x(1:end-1));
    y_dist = abs(y(2:end) - y(1:end-1));

    xy_dist = sqrt(x_dist.^2 + y_dist.^2);
end 

% some sort of smoothing, getting rid of extra points (assume less distance means less important) 
function [ret_x, ret_y] = resample_char(x,y,S,I)
    x = flip(x);
    y = flip(y);

    D = 0;
    p1 = [x(1),y(1)];
    ret_p = [p1];
    for i = 2:1:(size(x,1))
        p2 = [x(i),y(i)];
        dist = sqrt((p2(1) - p1(1))^2 + (p2(2) - p1(2))^2);
        num_seg = dist / S;
        if (num_seg < 1)
            continue;
        end
        if (num_seg == 1)
            ret_p = [ret_p; p2]; % append
            p1 = p2;
            continue;
        end
        num_seg = floor(num_seg);
        for j = 1:1:size(num_seg)
            dist = sqrt((p2(1) - p1(1))^2 + (p2(2) - p1(2))^2);
            % linear interpolate w p2
            p_eff(1) = p1(1) + S/dist * (p2(1) - p1(1));
            p_eff(2) = p1(2) + S/dist * (p2(2) - p1(2));
            ret_p = [ret_p; p_eff]; % append
            p1 = p_eff;
        end
    end % for each point in x
    ret_x = ret_p(:,1);
    ret_y = ret_p(:,2);
    
    ret_x = flip(ret_x);
    ret_y = flip(ret_y);
end

function [x_scaled,y_scaled] = scale_xy(x,y,a,b)
    x_scaled = (b - a) .* (x - min(x)) ./ (max(x) - min(x)) + a;
    y_scaled = (b - a) .* (y - min(y)) ./ (max(y) - min(y)) + a;
end

function R = rotx(degrees)
    a = degrees * pi/180;
    R = [1,0,0 ; 0,cos(a),-sin(a) ; 0,sin(a),cos(a)];
end

function R = rotz(degrees)
    a = degrees * pi/180;
    R = [cos(a),-sin(a),0 ; sin(a),cos(a),0 ; 0,0,1];
end

function R = roty(degrees)
    a = degrees * pi/180;
    R = [cos(a),0,sin(a) ; 0,1,0 ; -sin(a),0,cos(a)];
end
