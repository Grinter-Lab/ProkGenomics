while true; do
 wget -T 15 -c https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz && break
done

tar -xvzf gtdbtk_r214_data.tar.gz

