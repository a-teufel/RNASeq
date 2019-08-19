#! /bin/bash 
#! Prompts for README and runs main project script.
#! echo "Type README: "
#! read README

#! echo -n $README > README.txt
#sudo nohup /usr/local/MATLAB/R2016a/bin/glnxa64/MATLAB -singleCompThread -r "cd(pwd);addpath(genpath('../'));disp(pwd);ProjectScript"
#nohup /usr/bin/matlab -singleCompThread -r "cd(pwd);addpath(genpath('../'));disp(pwd);ProjectScript"
nohup /stor/system/opt/MATLAB/R2015b/bin/matlab -singleCompThread -r "cd(pwd);addpath(genpath('../'));disp(pwd);ProjectScript"
