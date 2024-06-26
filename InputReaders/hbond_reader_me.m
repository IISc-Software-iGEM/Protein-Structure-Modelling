%% 
% edited the file hbond_reader
%%
function hbond_reader_me(file_name,write_file_name)

if nargin ~= 2
    error('only two input arguments')
end

fid = fopen(file_name);
if fid < 0
    error('hbond file could not be opend')
end

%raw_data = textscan(fid,'%d %s %s %d %s %s');
raw_data = textscan(fid,'%d %s %s %d %s %s %f');
fclose(fid);

max_i = size(raw_data{1},1);
fid = fopen(write_file_name,'wt');
for i = 1:max_i
%     if strcmp(raw_data{3}{i},'O') && strcmp(raw_data{6}{i},'H')
%         upper_bound = 2;
%     elseif strcmp(raw_data{3}{i},'O') && strcmp(raw_data{6}{i},'N')
%         upper_bound = 3;
%    end
    
    if ( strcmp(raw_data{3}{i},'O') && strcmp(raw_data{6}{i},'H') ) || ...
            ( strcmp(raw_data{6}{i},'O') && strcmp(raw_data{3}{i},'H') )
        upper_bound = 2;
    elseif (strcmp(raw_data{3}{i},'O') && strcmp(raw_data{6}{i},'N')) || ...
              (strcmp(raw_data{6}{i},'O') && strcmp(raw_data{3}{i},'N'))
        upper_bound = 3;
    end    
    
    %fprintf(fid,'%4.0f  %s %s \t %4.0f %s %s \t %d \t \t #peak \t -1\n',...
    %    raw_data{1}(i),raw_data{2}{i},raw_data{3}{i},raw_data{4}(i),...
    %    raw_data{5}{i},raw_data{6}{i},upper_bound);
    
    fprintf(fid,'%4.0f  %s %s \t %4.0f %s %s \t %f \t \t #peak \t -1\n',...
        raw_data{1}(i),raw_data{2}{i},raw_data{3}{i},raw_data{4}(i),...
        raw_data{5}{i},raw_data{6}{i},raw_data{7}(i));
end
fclose(fid);

end

