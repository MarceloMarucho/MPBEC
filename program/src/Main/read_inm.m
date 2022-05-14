function [dime, glen, wout, bulk, bc,digits, selectsolver, dielx_str, diely_str, dielz_str, kappa_str, charge_str, pqr_str, pqr_cent_str, srfm_str, chgm_str, ener_str,in_nam_str, out_nam_str] = read_inm(inputfile)

file_id = fopen(inputfile, 'rt'); % open the file

%% dime
line=fgetl(file_id);
[dime, count, errmsg, nextindex]=sscanf(line, '%e');
dime=dime';
 
%% glen
line=fgetl(file_id);
[glen, count, errmsg, nextindex]=sscanf(line, '%e');
glen=glen';

%% bulk properties
line=fgetl(file_id);
[bulk, count, errmsg, nextindex]=sscanf(line, '%e');
bulk=bulk';

%% bc
bc=fgetl(file_id);

%% significant digits of rpecision
line=fgetl(file_id);
[strg, count, errmsg, nextindex]=sscanf(line, '%s%c');
[digits selectsolver]=strread(strg, '%d %s');
%% dielx_str
dielx_str=fgetl(file_id);

%% diely_str
diely_str=fgetl(file_id);

%% dielz_str
dielz_str=fgetl(file_id);

%% kappa_str
kappa_str=fgetl(file_id);

%% charge_str
charge_str=fgetl(file_id);

%% pqr_str
pqr_str=fgetl(file_id);

%% pqr_cent_str
pqr_cent_str=fgetl(file_id);

%% srfm_str
srfm_str=fgetl(file_id);

%% chgm_str
chgm_str=fgetl(file_id);

%% energy calculation
ener_str=fgetl(file_id);

%% input path nam_str
in_nam_str=fgetl(file_id);

%% output path nam_str
out_nam_str=fgetl(file_id);

%% writeout files
line=fgetl(file_id);
[wout, count, errmsg, nextindex]=sscanf(line, '%e');
wout=wout';

fclose(file_id); % close the file