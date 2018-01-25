#sourceme

# submit a bunch of jobs to scan xsec and mX using
# SLURM batch manager.

output_dir=fits-scan
batch_partition=hep
trials_per_job=50
timelimit=20

mkdir -p $output_dir

for i in {0..99}; do
	for xsec in 0.0 0.25 0.5 0.75 1.0; do
		for mX in {260..450..10}; do
			sbatch -o /dev/null -p $batch_partition -t $timelimit \
			./bias-test.py \
				--ws 3000_ggbb_lowmass.root \
				--freeze-bias \
				--xsec ${xsec} \
				--ntrial $trials_per_job \
				--seed $((mX+i)) \
				--mX $mX \
				--only-good \
				--out $output_dir/fits-x${xsec}-m${mX}-${i}.npy;
		done
	done
done
