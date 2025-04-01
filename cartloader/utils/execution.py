import os, sys, subprocess
from datetime import datetime
import random

# Execute makefile by SLURM using a command
def submit_cmd(mkpref, jobname, slurm, submit):
    slurm_account=slurm.get("account", None)
    slurm_partition=slurm.get("partition", None)
    slurm_mail_user=slurm.get("mail_user", None)
    slurm_cpus_per_task=slurm.get("cpus_per_task", None)
    slurm_mem=slurm.get("mem", None)
    slurm_hours=slurm.get("hours", None)
    if submit:
        submit_cmd = " ".join([
            "sbatch",
            f"--account={slurm_account}",
            f"--partition={slurm_partition}",
            f"--mail-user={slurm_mail_user}",
            "--mail-type=END,FAIL,REQUEUE",
            f"--job-name={jobname}",
            f"--output={mkpref}_%j.out",
            f"--error={mkpref}_%j.err",
            f"--time={slurm_hours}:00:00",
            "--nodes=1",
            "--ntasks=1",
            f"--cpus-per-task={slurm_cpus_per_task}",
            f"--mem={slurm_mem}",
            f"--wrap=\"make -f {mkpref}.mk -j 10\""
        ])
        print("Submitting job with command:", submit_cmd)
        subprocess.run(submit_cmd, shell=True)


# Generate a job file for local or SLURM
def write_jobfile(job_id, log_dir, cmds, slurm, submit_mode):
    jobpref = os.path.join(log_dir, job_id)
    # Generate SLURM script content
    if submit_mode == "local":
        job_content = """#!/bin/bash\n""" + "\n".join(cmds)
    else:
        # ?? The minimum required SLURM parameters
        missing_params = [k for k in [ "cpus_per_task", "mem_per_cpu"] if not getattr(slurm, k, None)]
        if missing_params:
            print(f"Error: Missing required SLURM parameters: {', '.join(missing_params)}")
            return
        # Generate SLURM options
        slurm_options = {
            "account": f"#SBATCH --account={slurm.account}" if getattr(slurm, "account", None) else "",
            "partition": f"#SBATCH --partition={slurm.partition}" if getattr(slurm, "partition", None) else "",
            "mail_user": f"#SBATCH --mail-user={slurm.mail_user}" if getattr(slurm, "mail_user", None) else "",
            "mail_type": "#SBATCH --mail-type=END,FAIL,REQUEUE" if getattr(slurm, "mail_user", None) else "",
            "job_name": f"#SBATCH --job-name={job_id}",
            "output": f"#SBATCH --output={jobpref}_%j.out",
            "error": f"#SBATCH --error={jobpref}_%j.err",
            "time": f"#SBATCH --time={getattr(slurm, 'hours', '5')}:00:00",
            "nodes": "#SBATCH --nodes=1",
            "ntasks": "#SBATCH --ntasks=1",
            "cpus_per_task": f"#SBATCH --cpus-per-task={getattr(slurm, 'cpus_per_task', 1)}",
            "mem_per_cpu": f"#SBATCH --mem-per-cpu={getattr(slurm, 'mem_per_cpu', '6500mb')}"
        }
        slurm_header = ["#!/bin/bash"] + [value for value in slurm_options.values() if value]  # Filter out empty entries and construct the header
        job_content = "\n".join(slurm_header + cmds)

    # Write SLURM script
    job_file = f"{jobpref}.job"
    try:
        with open(job_file, "w") as f:
            f.write(job_content)
        #print(f"A job file created: {job_file}")
        # run chmod 755 on the job file to make it executable
        os.chmod(job_file, 0o755) # 0o755
        return job_file
    except Exception as e:
        print(f"Error writing a job file: {e}")
        return
    

# Submit a job file to local or SLURM
def submit_job(job_file, submit_mode):
    try:
        if submit_mode == "local":
            result = subprocess.run(f"bash {job_file}", shell=True, check=True, capture_output=True, text=True)
            print("Job executed successfully:\n", result.stdout)
            # generate a random integer with 6 digits
            rand = random.randint(100000, 999999)
            log_pref = job_file.replace(".job", f"_{rand}")
            with open(f"{log_pref}.out", "w") as f:
                f.write(result.stdout)
            with open(f"{log_pref}.err", "w") as f:
                f.write(result.stderr)
        elif submit_mode == "slurm":
            result = subprocess.run(f"sbatch {job_file}", shell=True, check=True, capture_output=True, text=True)
            print("Job submitted successfully:\n", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error submitting job:\n", e.stderr)
