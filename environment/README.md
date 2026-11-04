# Runtime and resource planning

The recorded calculation used Python 3.11.13, PySCF 2.9.0, and block2-mpi 0.5.3. See `requirements.txt`; remaining dependencies are listed but not claimed to be a complete locked environment. The actual chemistry environment must match your cluster's MPI/BLAS runtime. Packaging and dry-run checks need only standard-library Python 3.11.

If your site provides a script that loads this environment, point `FE2S2_RUNTIME` at it:

```bash
export FE2S2_RUNTIME=/path/to/your/chemistry_env.sh
```

Otherwise activate your own compatible environment. `runtime.sh` optionally sources `FE2S2_RUNTIME`, then sets thread variables. No Python environment is silently installed or replaced.

Submit from the repository root. The batch wrapper uses `SLURM_SUBMIT_DIR`, not its temporary Slurm spool location:

```bash
sbatch --cpus-per-task=8 --mem=64G --time=02:00:00   environment/submit_step.sh 02_parent_uks --workdir /scratch/fe2s2_example

sbatch --cpus-per-task=16 --mem=160G --time=12:00:00   environment/submit_step.sh 09_production_states --workdir /scratch/fe2s2_example   --action refine-cas20 --bonds 3000 4000 5000
```

Initialize the work directory with step 01 first. The template sets no partition or account; add `--partition` (and `--account` if your site requires it) to the `sbatch` command. Template defaults are 16 CPUs, 160 GiB, 12 hours; use each step's `step.json`/README to avoid over-requesting small stages. The wrapper submits no dependent jobs automatically. Inspect completion and convergence before submitting the next step.

The DMRG driver uses one process with threaded block2 (`mpi=False`), although the installed distribution is block2-mpi. Steps 04/05 default to a 32 GiB block2 stack; step 09 defaults to 80 GiB. `DMRG_STACK_GB` can override these, within your allocation. FCI/direct validation should have at least the documented 64 GiB allocation; `VALIDATION_MEMORY_MB` defaults to 24000 for step 10.

Budget at least 50 GiB for durable production arrays alone and substantial additional space for DMRG checkpoints, MPO scratch, and refinement backups. A practical scratch budget of several hundred GiB is prudent. The full factor file is ~26.23 GiB and AO packed factors ~15.07 GiB. Do not materialize a full retained four-index ERI tensor (~1.06 TiB) or the dense AO tensor (~1.40 TiB).

For validation of an installed chemistry environment, run step 10's `smoke` action in a small compute allocation after step 01. No molecular recalculation was performed merely to package this repository.
