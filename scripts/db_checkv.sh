rm -r checkv-db*

while true;do
    wget -T 15 -c https://portal.nersc.gov/CheckV/checkv-db-v0.6-full.tar.gz && break
done

tar -zxvf checkv-db-v0.6-full.tar.gz  
cd checkv-db-v0.6/genome_db
diamond makedb --in checkv_reps.faa --db checkv_reps