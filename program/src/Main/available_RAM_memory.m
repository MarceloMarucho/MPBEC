function [ availmemm ] = available_RAM_memory() 
% which operating system?

strp=computer;

 if (strcmp(strp,'PCWIN')==1 ||strcmp(strp,'PCWIN64')==1)
% read available memory from terminal and export data in text file
system('systeminfo | find "MB" > usedmemory.txt');
pause(1);
% read data from file
file_id = fopen('usedmemory.txt', 'rt'); % open the file
line1=fgetl(file_id);
line=fgetl(file_id);
% delete special (non-alphanumeric) characters
expression = ',';
replace = '';
TT = regexprep(line, expression, replace);
if (strcmp(TT,'')==1)
    TT=line;
end
expression = '.';
replace = '';
T = strrep(TT, expression, replace);
if (strcmp(T,'')==1)
    T=TT;
end
str = {char(128), char(131), char(138), char(142), char(154), ...
      char(158), char(159), char(162), char(165), char(181), ...
      char(192), char(193), char(194), char(195), char(196), ...
      char(197), char(199), char(200), char(201), char(202), ...
      char(203), char(204), char(205), char(206), char(207), ...
      char(209), char(210), char(211), char(212), char(213), ...
      char(214), char(216), char(217), char(218), char(219), ...
      char(220), char(221), char(224), char(225), char(226), ...
      char(227), char(228), char(229), char(231), char(232), ...
      char(233), char(234), char(235), char(236), char(237), ...
      char(238), char(239), char(241), char(242), char(243), ...
      char(244), char(245), char(246), char(248), char(249), ...
      char(250), char(251), char(252), char(253), char(255),char(161)};
strreplace = {'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U','o'};
    for i = 1:numel(str)
        T = regexprep(T, str{i}, strreplace{i});
    end
% get only the number of available memory (in Mbs)
memm2 = strread(lower(T), '%s', 'whitespace',['a':'z' ':' ' \t']);
fclose(file_id); 
delete('usedmemory.txt')
% checking the data
if (numel(memm2)~= 1)
disp('Sorry. MPBEC was unable to determine the available RAM memory. This may be due to operating system language issues.')
    disp(' ')
    availmemm=-1;
    return
elseif ischar(str2num(cell2mat(memm2(1))))
disp('Sorry. MPBEC was unable to determine the available RAM memory. This may be due to operating system language issues.')
    disp(' ')
    availmemm=-1;
    return
end
memm=str2num(cell2mat(memm2(1)));
availmemm=memm;
% alert message if available RAM memory is lower than 100Mbs.
if (availmemm<100.0)
 argg2=sprintf('%s %g %s','the available RAM memory in your computer is low (',availmemm,' Mbs). Please consider to abort the calculations, close unsed programs, free cached RAM memory and try it again. If the problem persists, you would need to increase the RAM memory of your computer. Do you want to abort the calculations?? ');
button= questdlg(argg2, ...
        'Exit Dialog','Yes','No','Yes');
         switch button
            case 'Yes',
                disp(' ');
              disp('Exiting MPBEC calculation. This may take a few seconds...');
              exitmpbec = 1;
            case 'No',
                disp(' ');
               disp('resume calculations...');
         end
end
 end

% similar explanation for Linux and Mac

if (strcmp(strp,'GLNX86')==1 ||strcmp(strp,'GLNXA64')==1)
system('free  -m | grep Mem: > usedmemory.txt');
pause(1);
file_id = fopen('usedmemory.txt', 'rt'); 
line=fgetl(file_id);

expression = ',';
replace = '';
TT = regexprep(line, expression, replace);
if (strcmp(TT,'')==1)
    TT=line;
end
expression = '.';
replace = '';
T = strrep(TT, expression, replace);
if (strcmp(T,'')==1)
    T=TT;
end
str = {char(128), char(131), char(138), char(142), char(154), ...
      char(158), char(159), char(162), char(165), char(181), ...
      char(192), char(193), char(194), char(195), char(196), ...
      char(197), char(199), char(200), char(201), char(202), ...
      char(203), char(204), char(205), char(206), char(207), ...
      char(209), char(210), char(211), char(212), char(213), ...
      char(214), char(216), char(217), char(218), char(219), ...
      char(220), char(221), char(224), char(225), char(226), ...
      char(227), char(228), char(229), char(231), char(232), ...
      char(233), char(234), char(235), char(236), char(237), ...
      char(238), char(239), char(241), char(242), char(243), ...
      char(244), char(245), char(246), char(248), char(249), ...
      char(250), char(251), char(252), char(253), char(255),char(161)};
strreplace = {'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U','o'};
    for i = 1:numel(str)
        T = regexprep(T, str{i}, strreplace{i});
    end
memm = strread(lower(T), '%s', 'whitespace',['a':'z' ':' ' \t']);
fclose(file_id); 
delete('usedmemory.txt')
if (numel(memm)~= 6)
    disp('Sorry. MPBEC was unable to determine the available RAM memory. This may be due to operating system language issues.')
    disp(' ')
availmemm=-1;
    return
elseif ischar(str2num(cell2mat(memm(3))))
    disp('Sorry. MPBEC was unable to determine the available RAM memory. This may be due to operating system language issues .')
    disp(' ')
availmemm=-1;
    return
end

availmemm=str2num(cell2mat(memm(3)));

if (availmemm<100.0)
 argg2=sprintf('%s %g %s','the available RAM memory in your computer is low (',availmemm,' Mbs). Please consider to abort the calculations, close unsed programs, free cached RAM memory and try it again. If the problem persists, you would need to increase the RAM memory of your computer. Do you want to abort the calculations?? ');
button= questdlg(argg2, ...
        'Exit Dialog','Yes','No','Yes');
         switch button
            case 'Yes',
                disp(' ');
              disp('Exiting MPBEC calculation. This may take a few seconds...');
              exitmpbec = 1;
            case 'No',
                disp(' ');
               disp('resume calculations...');
         end
end
end


if (strcmp(strp,'MACI')==1 ||strcmp(strp,'MACI64')==1)
system('top -l 1 | head -n 10 | grep PhysMem > usedmemory.txt');
pause(1);
file_id = fopen('usedmemory.txt', 'rt'); 

line=fgetl(file_id);

expression = ',';
replace = '';
TT = regexprep(line, expression, replace);
if (strcmp(TT,'')==1)
    TT=line;
end
expression = '.';
replace = '';
T = strrep(TT, expression, replace);
if (strcmp(T,'')==1)
    T=TT;
end
str = {char(128), char(131), char(138), char(142), char(154), ...
      char(158), char(159), char(162), char(165), char(181), ...
      char(192), char(193), char(194), char(195), char(196), ...
      char(197), char(199), char(200), char(201), char(202), ...
      char(203), char(204), char(205), char(206), char(207), ...
      char(209), char(210), char(211), char(212), char(213), ...
      char(214), char(216), char(217), char(218), char(219), ...
      char(220), char(221), char(224), char(225), char(226), ...
      char(227), char(228), char(229), char(231), char(232), ...
      char(233), char(234), char(235), char(236), char(237), ...
      char(238), char(239), char(241), char(242), char(243), ...
      char(244), char(245), char(246), char(248), char(249), ...
      char(250), char(251), char(252), char(253), char(255),char(161)};
strreplace = {'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U', ...
    'A', 'E', 'I', 'O', 'U','o'};
    for i = 1:numel(str)
        T = regexprep(T, str{i}, strreplace{i});
    end
memm = strread(lower(T), '%s', 'whitespace',['a':'z' ':' '(' ')' ',' ' \t']);
fclose(file_id); 
delete('usedmemory.txt')
if (numel(memm)~= 3)
disp('Sorry. MPBEC was unable to determine the available RAM memory. This may be due to operating system language issues.')
    disp(' ')
    availmemm=-1;
    return
elseif ischar(str2num(cell2mat(memm(3))))
disp('Sorry. MPBEC was unable to determine the available RAM memory. This may be due to operating system language issues.')
    disp(' ')
    availmemm=-1;
    return
end

availmemm=str2num(cell2mat(memm(3)));
if (availmemm<100.0)
 argg2=sprintf('%s %g %s','the available RAM memory in your computer is low (',availmemm,' Mbs). Please consider to abort the calculations, close unsed programs, free cached RAM memory and try it again. If the problem persists, you would need to increase the RAM memory of your computer. Do you want to abort the calculations?? ');
button= questdlg(argg2, ...
        'Exit Dialog','Yes','No','Yes');
         switch button
            case 'Yes',
                disp(' ');
              disp('Exiting MPBEC calculation. This may take a few seconds...');
              exitmpbec = 1;
            case 'No',
                disp(' ');
               disp('resume calculations...');
         end
end
end

clear strp memm



