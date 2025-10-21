mkdir -p profiles
python get_profile.py data/e_coli_HiSeq2K/e_coli_CNR0028307_p33_1.fastq.gz profiles/pe100_real_1.profile 100 &
python get_profile.py data/e_coli_HiSeq2K/e_coli_CNR0028307_p33_2.fastq.gz profiles/pe100_real_2.profile 100 &

python get_profile.py data/soybean_HiSeq2500/soybean_SRR16074289_1.fastq.gz profiles/pe300_real_1.profile 300 &
python get_profile.py data/soybean_HiSeq2500/soybean_SRR16074289_2.fastq.gz profiles/pe300_real_2.profile 300 &
wait
