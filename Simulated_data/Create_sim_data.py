import subprocess
import os
import random

# ================= User Configuration =================
# Ensure the reference genome index (.fai) exists in the directory
REF_GENOME = "Homo_sapiens.GRCh37.dna.primary_assembly.fa"
OUTPUT_DIR = "."
DEPTH = 40                   # Sequencing Depth for fusion regions
READ_LEN = 100               # Read Length
BACKGROUND_READS = 100000    # Number of background reads
# ======================================================

def run_cmd(cmd):
    # print(f"[Run] {cmd}")
    subprocess.check_call(cmd, shell=True)

def introduce_snps(sequence, rate=0.05):
    """
    Introduce variants for Group C (High SNP).
    Rate is set to 0.05 (5%) to challenge alignment tools but allow assembly.
    """
    seq = list(sequence)
    # Only introduce mutations in the middle 40% (near the breakpoint)
    start = int(len(seq) * 0.3)
    end = int(len(seq) * 0.7)
    num_muts = int((end - start) * rate)
    indices = random.sample(range(start, end), num_muts)
    for i in indices:
        choices = [b for b in "ACGT" if b != seq[i]]
        seq[i] = random.choice(choices)
    return "".join(seq)

def generate_dataset():
    # === The Corrected Essential 9 List (GRCh37/Ensembl Style - No 'chr') ===
    fusions = [
        # --- Group A: Standard (AF 0.5) ---
        {"name": "Std_EWSR1_FLI1",  "chrA": "22", "posA": 29683123,   "chrB": "11", "posB": 128573500, "af": 0.5, "type": "clean"},
        {"name": "Std_BCR_ABL1",    "chrA": "22", "posA": 23524426,   "chrB": "9",  "posB": 133729451, "af": 0.5, "type": "clean"},
        {"name": "Std_TMPRSS2_ERG", "chrA": "21", "posA": 42870000,   "chrB": "21", "posB": 39775000,  "af": 0.5, "type": "clean"},

        # --- Group B: Low VAF (AF 0.2) ---
        {"name": "Low_EML4_ALK",    "chrA": "2",  "posA": 42522656,   "chrB": "2",  "posB": 29446394,  "af": 0.2, "type": "clean"},
        {"name": "Low_CD74_ROS1",   "chrA": "5",  "posA": 149782000,  "chrB": "6",  "posB": 117642000, "af": 0.2, "type": "clean"},
        {"name": "Low_KMT2A_AFF1",  "chrA": "11", "posA": 118353000,  "chrB": "4",  "posB": 88000000,  "af": 0.2, "type": "clean"},

        # --- Group C: High SNP (AF 0.5) ---
        {"name": "SNP_FGFR3_TACC3", "chrA": "4",  "posA": 1806000,    "chrB": "4",  "posB": 1737000,   "af": 0.5, "type": "snp"},
        {"name": "SNP_RUNX1_RUNX1T1","chrA": "21","posA": 36252000,   "chrB": "8",  "posB": 93058000,  "af": 0.5, "type": "snp"},
        {"name": "SNP_NPM1_ALK",    "chrA": "5",  "posA": 170837543,  "chrB": "2",  "posB": 29443000,  "af": 0.5, "type": "snp"},
    ]

    r1_list, r2_list = [], []
    print(f">>> Start generating data: ReadLen={READ_LEN}bp, Depth={DEPTH}x")

    # 1. Generate Fusions
    for f in fusions:
        print(f"   Processing: {f['name']}")
        # Extract sequences (600bp flanking)
        # Note: 'chr' prefix removed from chromosome names
        run_cmd(f"samtools faidx {REF_GENOME} {f['chrA']}:{f['posA']-600}-{f['posA']} > tA.fa")
        run_cmd(f"samtools faidx {REF_GENOME} {f['chrB']}:{f['posB']}-{f['posB']+600} > tB.fa")
        
        with open(f"{OUTPUT_DIR}/{f['name']}.fa", 'w') as out:
            out.write(f">{f['name']}\n")
            with open("tA.fa") as fA: sA = "".join(fA.readlines()[1:]).strip()
            with open("tB.fa") as fB: sB = "".join(fB.readlines()[1:]).strip()
            final = sA + sB
            if f['type'] == 'snp': 
                # Introduce SNPs at 5% rate for Group C
                final = introduce_snps(final, rate=0.05) 
            out.write(final + "\n")
        
        cov = max(2, int(DEPTH * f['af']))
        pfx = f"{OUTPUT_DIR}/{f['name']}"
        # ART params: -l 100 (Read Length), -f (Coverage)
        run_cmd(f"art_illumina -ss HS25 -p -l {READ_LEN} -f {cov} -m 350 -s 50 -i {OUTPUT_DIR}/{f['name']}.fa -o {pfx} -q")
        r1_list.append(f"{pfx}1.fq"); r2_list.append(f"{pfx}2.fq")

    # 2. Generate Background 
    print(f">>> Generating background reads ({BACKGROUND_READS})...")
    # Using Chromosome 3 (Ensembl style '3') for background
    bg_start = 50000000
    run_cmd(f"samtools faidx {REF_GENOME} 3:{bg_start}-{bg_start+2000000} > {OUTPUT_DIR}/bg_ref.fa")
    
    bg_pfx = f"{OUTPUT_DIR}/background"
    # ART params: -c (Count, specific number of reads)
    run_cmd(f"art_illumina -ss HS25 -p -l {READ_LEN} -c {BACKGROUND_READS} -m 350 -s 50 -i {OUTPUT_DIR}/bg_ref.fa -o {bg_pfx} -q")
    r1_list.append(f"{bg_pfx}1.fq"); r2_list.append(f"{bg_pfx}2.fq")

    # 3. Merge
    print(">>> Merging files...")
    run_cmd(f"cat {' '.join(r1_list)} > test_R1.fq")
    run_cmd(f"cat {' '.join(r2_list)} > test_R2.fq")
    
    # Cleanup
    # Note: Removed the 'rm' command for safety to avoid deleting reference genome.
    # Please manually remove the temporary .fa and .fq files (excluding test_R1/R2.fq).
    
    print(f"\n✅ Done! Dataset generated:\n - test_R1.fq\n - test_R2.fq\n (Read Length: {READ_LEN}, Background: {BACKGROUND_READS} reads)")

if __name__ == "__main__":
    generate_dataset()
