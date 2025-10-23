# Benchmarks

Benchmark of `art_modern` with WGSim (Bundled with SAMtools) and the original ART simulator. All benchmarks are performed on a Linux machine with Intel DPC++/C++ compiler (version `2025.2.1 (2025.2.0.20250806)`) with `-O3` and `-mtune=native`.

The GCC build of `art_modern` (Using GCC version `13.3.0-6ubuntu2~24.04`) is also included for reference.

`art_modern` was linked to Jemalloc version `5.3.0-0-g54eaed1d8b56b1aa528be3bdd1877e59c56fa90c`, MKL version `2025.0.2.0 (20250200)` and Boost version `1.83.0`.

## Instructions

Use `build-others.sh` to build other simulators.

Use `build-update.sh` to build and update `art_modern`.

Use `run-*.sh` to run benchmark.

## Datasets

1. Reference genome of _C. Elegans_ with 10X coverage.
2. Long mRNAs of humans that exceed 500 nucleotides with 4X coverage.
3. Large contig dataset: Simulated genome with 4 contigs, each containing 32 G bases.
4. Large mRNAs dataset: Simulated genome with 32 G contigs, each containing 500 K bases.

The simulators were asked to generate a FASTQ on both datasets with paired-end reads lengthen 100bp.

## Results

### PC 1

- **MODEL** HP ZBook Power 15.6 inch G10 Mobile Workstation PC [spec](https://support.hp.com/us-en/document/ish_8090133-8090177-16)
- **CPU** 13th Gen Intel(R) Core(TM) i7-13700H [spec](https://www.intel.com/content/www/us/en/products/sku/232128/intel-core-i713700h-processor-24m-cache-up-to-5-00-ghz/specifications.html)
- **MEM** 2 * 16GB Hynix/Hyundai 5600MHz No-ECC (HMCG78AGBSA095N)
- **SSD** WD_BLACK SN750 NVMe SSD (WDS200T3X0C-00SJG0)
- **OS** Linux Mint 22
- **KERNEL** Linux 6.8.0-79-generic

## Used Third-Party Codes

### SAMtools 1.21 Release Tarball

Available from <https://sourceforge.net/projects/samtools/files/samtools/1.21/samtools-1.21.tar.bz2>

Affected files:

- `src/wgsim.c`

With MIT/Expat License.

### Original ART Source Code v2016.06.05

Available from <https://www.niehs.nih.gov/sites/default/files/2024-02/artsrcmountrainier2016.06.05linux.tgz>

Affected files:

- `src/art_original/*`

With GNU GPL v3 License.

### DWGSIM 0.1.15 Repo

Available from <https://github.com/nh13/DWGSIM>, commit `163fede`.

Affected files:

- `src/dwgsim/*`

With GNU GPL v2 License.

### pIRS 2 Repo

Available from <https://git.code.sf.net/p/pirsim/code>, commit `bee9b594`.

Affected files:

- `src/pirs/*`

With GNU GPL v2 License. Some minor changes were made to comfort Intel compilers.
