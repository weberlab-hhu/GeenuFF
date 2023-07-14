mkdir three_algae
cd three_algae

# download all the data
for pth in ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/chlamydomonas_reinhardtii/dna/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.dna_sm.toplevel.fa.gz \
           ftp://ftp.ensemblgenomes.org/pub/plants/release-47/gff3/chlamydomonas_reinhardtii/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.47.gff3.gz \
ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/ostreococcus_lucimarinus/dna/Ostreococcus_lucimarinus.ASM9206v1.dna.toplevel.fa.gz \
           ftp://ftp.ensemblgenomes.org/pub/plants/release-47/gff3/ostreococcus_lucimarinus/Ostreococcus_lucimarinus.ASM9206v1.47.gff3.gz \
           ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/cyanidioschyzon_merolae/dna/Cyanidioschyzon_merolae.ASM9120v1.dna.toplevel.fa.gz \
           ftp://ftp.ensemblgenomes.org/pub/plants/release-47/gff3/cyanidioschyzon_merolae/Cyanidioschyzon_merolae.ASM9120v1.47.gff3.gz
do
  wget $pth
  sleep 0.4s
done

# uncompress the data
gunzip *.gz

# put the data in the format compatible with the --basedir parameter
# basically --basedir needs a folder <your_species>
# with a subfolder <your_species>/input containing a (compressed) gff3 annotation and fasta genome file.
# The results will then be located in <your_species>/output
# If desired, you can alternatively specify all file parameters individually
# --gff3 <your.gff3> --fasta <your.fa> --db-path <your_output_genuff.sqlite3> --log-file <your_output.log>
species="Chlamydomonas_reinhardtii Ostreococcus_lucimarinus Cyanidioschyzon_merolae"

for sp in $species
do
  spdir=$sp/input
  mkdir -p $spdir
  mv ${sp}.* $spdir/
done

# import into databases (the main output will land in <basedir>/output/<species>.sqlite3
for sp in $species
do
  import2geenuff.py --basedir $sp --species $sp
done
cd ..
