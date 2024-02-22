scp bio1:/home/lper0012/tasks/rhys.grinter/pipeline/main.nf ./
scp bio1:/home/lper0012/tasks/rhys.grinter/pipeline/modules ./
scp bio1:/home/lper0012/tasks/rhys.grinter/pipeline/*config* ./
git add * --all
git commit -m $1
git push origin main

