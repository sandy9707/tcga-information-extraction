# #!/bin/zsh
# cd '/Users/sandy/ResilioSync/git/YeziPublish'
# git pull
# git add .
# commitTime=`date +"%Y-%m-%d %H-%M-%S"`
# git commit -a -m "${commitTime}"
# git push origin main

SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
echo $SHELL_FOLDER
cd $SHELL_FOLDER
git add .
commitTime=`date +"%Y-%m-%d %H-%M-%S"`
git commit -a -m "${commitTime}"
git push origin main