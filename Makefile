
run:
	cargo run -- to-fastq \
	--i1 /Users/mahogny/Downloads/joram/Joram-singlecell-2nd_S1_L001_R1_001.fastq.gz \
	--i2 /Users/mahogny/Downloads/joram/Joram-singlecell-2nd_S1_L001_R2_001.fastq.gz \
	--o1 /Users/mahogny/Downloads/joram/rustout_R1.fastq.gz \
	--o2 /Users/mahogny/Downloads/joram/rustout_R2.fastq.gz \
	--h /Users/mahogny/Downloads/joram/rustout_hist.csv

help:
	cargo run -- -help

getcmake:
	/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
	brew install cmake
