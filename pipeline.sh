



reference=$1
assembly_contigs=$2
scaffolds=$3
outdir=$4


mkdir -p $outdir


#echo $assembly_contigs $reference "$outdir/ref_scaf"


# PRODUCE PERFECT CONTIGS

#python3 build_ref_scaf.py assembly_contigs.fa reference.fa ref_scaf
python3 build_ref_scaf.py $assembly_contigs $reference "$outdir/ref_scaf"


python3 build_out_scaf.py $scaffolds "$outdir/ref_scaf.fa" "$outdir/out_scaf"
#python3 build_out_scaf.py scaffmatch.fa ref_scaf.fa out_scaf


# RUN VALIDATION

python validation.py $outdir/out_scaf.scaf $outdir/ref_scaf.scaf




