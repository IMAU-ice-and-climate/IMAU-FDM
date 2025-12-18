filename = 'T_firn_DYE-2_16.nc';
finfo = ncinfo(filename);

disp('Displaying variables')
disp({finfo.Variables.Name})
disp('Displaying attributes linked to time variables')
tmp = {finfo.Variables.Attributes};
disp({tmp{1}.Name})
disp({tmp{1}.Value})

% extracting
time = ncread(filename,'time'); %time in  'days since 1900-1-1 0:0:0' 
% changing to matlab's native time reference (in days since year 0 day 0)
time_matlab = time + datenum(1900,1,1);

disp('Displaying first three obsevation times and last one')
disp(datestr(time_matlab([1:3 end])))
