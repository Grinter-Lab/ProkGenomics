rm -r checkv-db*
curl -C --retry 5 --retry-all-errors https://portal.nersc.gov/CheckV/checkv-db-v0.6-full.tar.gz
tar -zxvf checkv-db-v0.6-full.tar.gz  
cd checkv-db-v0.6
diamond makedb --in checkv_reps.faa --db checkv_reps