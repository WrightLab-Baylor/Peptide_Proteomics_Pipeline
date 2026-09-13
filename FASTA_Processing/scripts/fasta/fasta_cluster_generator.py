#!/usr/bin/env python3
import argparse
import subprocess
import os
import shutil
from Bio import SeqIO
from glob import glob
from pathlib import Path
import sys
import re
import tempfile

# Resolve cd-hit from --cdhit-path or PATH only.

def str2bool(val):
    return str(val).lower() in ("yes", "true", "t", "1")

def find_default_input():
    fasta_files = glob("database/*.fasta") + glob("database/*.faa")
    if len(fasta_files) == 0:
        raise FileNotFoundError("❌ No FASTA/FAA file found in ./database/")
    elif len(fasta_files) > 1:
        raise RuntimeError("❌ Multiple FASTA/FAA files found in ./database/. Please use -i to specify input.")
    return fasta_files[0]

def is_already_processed(fasta_path, max_checks=1000):
    """Check up to `max_checks` FASTA entries for 'Cluster=' in header."""
    count = 0
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
                    if 'Cluster=' in line:
                        return True
                    if count >= max_checks:
                        break
    except Exception as e:
        print(f"[ERROR] Failed to check for existing cluster annotations: {e}")
    return False

def resolve_cdhit_path(provided_path=None):
    """
    Resolve cd-hit executable path.
    - If provided_path is set, validate it is an executable and return it.
    - Otherwise try to find 'cd-hit' on PATH via shutil.which.
    - Exit with an informative message if not found.
    """
    if provided_path:
        if os.path.isfile(provided_path) and os.access(provided_path, os.X_OK):
            return provided_path
        else:
            sys.exit(f"[ERROR] Provided cd-hit path '{provided_path}' is not an executable.")
    found = shutil.which("cd-hit")
    if found:
        return found
    sys.exit("[ERROR] cd-hit binary not found in PATH and no --cdhit-path provided. Install cd-hit or pass --cdhit-path.")

def run_cdhit_clustering(input_fasta, output_fasta, identity, word_len, threads, memory, cdhit_path):
    out_dir = os.path.dirname(os.path.abspath(output_fasta)) or "."
    os.makedirs(out_dir, exist_ok=True)

    cmd = [
        cdhit_path,
        "-i", input_fasta,
        "-o", output_fasta,
        "-c", str(identity),
        "-n", str(word_len),
        "-T", str(threads),
        "-M", str(memory),
        "-d", "0" #Prevents FASTA header truncation
    ]

    print(f"[INFO] Running cd-hit: {' '.join(cmd)}")
    # Stream output directly to terminal
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)
    
    #Confirm outputs exist
    clstr_path = output_fasta + ".clstr"
    if not os.path.exists(output_fasta) or not os.path.exists(clstr_path):
        raise FileNotFoundError(f"Expected cd-hit outputs not found: {output_fasta} or {clstr_path}")
    return Path(output_fasta), Path(clstr_path)

def parse_clstr_file(clstr_path):
    """
    Parse cd-hit .clstr file and return a mapping: sequence_id -> Cluster_XXXX
    """
    cluster_map = {}
    cluster_counter = 0
    current_cluster = None

    if not clstr_path.exists():
        raise FileNotFoundError(f"CD-HIT cluster file not found: {clstr_path}")

    with open(clstr_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                cluster_counter += 1
                current_cluster = f"Cluster_{cluster_counter:04d}"
                continue
            # member lines look like: "0    123aa, >seq_id... *"
            if current_cluster and '>' in line and '...' in line:
                part = line.split('>', 1)[1]
                seq_id = part.split('...')[0].strip()
                cluster_map[seq_id] = current_cluster
    return cluster_map

def candidate_seq_ids_for_record(rec):
    """
    Return a prioritized list of candidate sequence IDs to try when looking up cluster_map.
    Tries:
      1) accession if header is pipe-delimited (parts[1])
      2) rec.id
      3) rec.description (full header)
      4) variants with stripped version suffixes (e.g., /1) or spaces removed
    """
    candidates = []
    seq_id = rec.id
    # 1) accession-style (like >db|accession|rest)
    if '|' in seq_id:
        parts = seq_id.split('|')
        if len(parts) >= 2 and parts[1]:
            candidates.append(parts[1])
    # 2) raw id
    candidates.append(seq_id)
    # 3) full description (cd-hit may have used full header)
    if rec.description and rec.description != rec.id:
        candidates.append(rec.description)
    # 4) normalized variants
    # strip version suffixes like '/1' or '.1' at end
    def strip_version(s):
        for sep in ('/', '.'):
            if s.endswith(sep + "1") or s.endswith(sep + "2") or s.endswith(sep + "3"):
                return s.rsplit(sep, 1)[0]
        return s
    for cand in list(candidates):
        stripped = strip_version(cand)
        if stripped != cand:
            candidates.append(stripped)
        # also try removing whitespace
        ws_removed = cand.replace(" ", "")
        if ws_removed != cand:
            candidates.append(ws_removed)
    # dedupe while keeping order
    seen = set()
    ordered = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            ordered.append(c)
    return ordered

def _norm_uniprot(description: str) -> str | None:
    """
    Normalize UniProt headers like:
    sp|Contaminant_AAAA1|Anti-FLAG_AffinityTag ...
    Returns: 'Contaminant_AAAA1'
    """
    if not description:
        return None
    # Strip whitespace
    description = description.strip()
    # Extract only the first "word" (until first space) because id is in the first token
    first_token = description.split(" ")[0]
    # Case-insensitive match
    match = re.match(r"^[a-z]{2}\|([A-Z0-9_]+)\|", first_token, re.IGNORECASE)
    return match.group(1) if match else None

def _norm_pipe_accession(s):
    if '|' in s:
        parts = s.split('|')
        if len(parts) >= 2 and parts[1]:
            return parts[1]
    return None

def _norm_raw_id(s):
    return s

def _norm_description(s):
    return s  # caller will pass rec.description explicitly

def _norm_strip_version(s):
    # strip trailing /N or .N where N is small int
    for sep in ('/', '.'):
        if sep in s:
            base, tail = s.rsplit(sep, 1)
            if tail.isdigit():
                return base
    return s

def _norm_no_whitespace(s):
    return s.replace(" ", "")

# ordered list of named normalizers. Adjust priority if needed.
NORMALIZERS = [
    ("uniprot", lambda rec: _norm_uniprot(rec.description)),
    ("pipe_accession", lambda rec: _norm_pipe_accession(rec.id)),
    ("raw_id", lambda rec: _norm_raw_id(rec.id)),
    ("description", lambda rec: _norm_description(rec.description)),
    ("strip_version", lambda rec: _norm_strip_version(rec.id)),
    ("no_whitespace", lambda rec: _norm_no_whitespace(rec.id)),
]

def annotate_fasta_with_clusters(input_fasta, cluster_map, output_fasta, keep_temp=False):
    """
    Adaptive, streaming annotation:
      - Try a first-pass normalization determined from the first header.
      - Stream FASTA, write matched annotated records to output.
      - Unmatched records go into a temp FASTA.
      - Iteratively determine a new normalization from first unmatched, reprocess unmatched file.
      - Stop when no more unmatched can be matched or when no normalizer yields matches.
    keep_temp: if True, don't delete temporary unmatched files (for debugging).
    """
    def choose_normalizer_for_record(rec, cluster_keys):
        """Return normalizer function name and callable that yields a key if it matches cluster_keys, else None."""
        for name, func in NORMALIZERS:
            try:
                key = func(rec)
            except Exception:
                key = None
            if key and key in cluster_keys:
                return name, func
        return None, None

    cluster_keys = set(cluster_map.keys())

    # 1) sample first record to pick initial normalizer
    with open(input_fasta, "r") as fh:
        first_record = next(SeqIO.parse(fh, "fasta"), None)

    if first_record is None:
        raise ValueError("Input FASTA appears empty.")

    name, func = choose_normalizer_for_record(first_record, cluster_keys)
    if name is None:
        # fallback to raw_id as initial assumption
        name, func = "raw_id", lambda rec: rec.id

    print(f"[INFO] Initial normalizer chosen: {name}")

    # Prepare output file and initial temp unmatched file
    out_handle = open(output_fasta, "w")
    temp_in_path = input_fasta
    current_normalizer = func
    iteration = 0

    while True:
        iteration += 1
        # Process input FASTA (or previous unmatched) and write matches to output, unmatched to new temp file
        tmp_fd, temp_out_path = tempfile.mkstemp(suffix=".fasta", prefix="unmatched_")
        os.close(tmp_fd)
        unmatched_count = 0

        with open(temp_in_path, "r") as in_fh, open(temp_out_path, "w") as unmatched_fh:
            for rec in SeqIO.parse(in_fh, "fasta"):
                try:
                    key = current_normalizer(rec)
                except Exception:
                    key = None
                cluster_label = None
                if key and key in cluster_map:
                    cluster_label = cluster_map[key]
                else:
                    # fallback: also try rec.id and rec.description quickly (rare)
                    if rec.id in cluster_map:
                        cluster_label = cluster_map[rec.id]
                    elif rec.description in cluster_map:
                        cluster_label = cluster_map[rec.description]

                if cluster_label:
                    original_desc = rec.description
                    rec.description = f"{original_desc} Cluster={cluster_label.replace('Cluster_','')}"
                    SeqIO.write(rec, out_handle, "fasta")
                else:
                    SeqIO.write(rec, unmatched_fh, "fasta")
                    unmatched_count += 1

        print(f"[INFO] Iteration {iteration}: unmatched_count = {unmatched_count}")

        if unmatched_count == 0:
            # No more unmatched sequences, done
            break

        # Read first unmatched record to choose new normalizer
        with open(temp_out_path, "r") as ufh:
            first_unmatched = next(SeqIO.parse(ufh, "fasta"), None)

        if first_unmatched is None:
            # No unmatched records (should not happen here)
            break

        new_name, new_func = choose_normalizer_for_record(first_unmatched, cluster_keys)

        if new_name is None:
            # No suitable normalizer found, mark all remaining unmatched as Unassigned
            print(f"[INFO] No normalizer matched first unmatched (id={first_unmatched.id}). Marking remaining {unmatched_count} as Unassigned.")
            with open(temp_out_path, "r") as ufh:
                for rec in SeqIO.parse(ufh, "fasta"):
                    original_desc = rec.description
                    rec.description = f"{original_desc} Cluster=Unassigned"
                    SeqIO.write(rec, out_handle, "fasta")
            break

        if new_func == current_normalizer:
            # No progress with a new normalizer, avoid infinite loop
            print("[INFO] New normalizer equals current; no progress. Marking remaining unmatched as Unassigned.")
            with open(temp_out_path, "r") as ufh:
                for rec in SeqIO.parse(ufh, "fasta"):
                    original_desc = rec.description
                    rec.description = f"{original_desc} Cluster=Unassigned"
                    SeqIO.write(rec, out_handle, "fasta")
            break

        # Prepare for next iteration
        print(f"[INFO] Switching normalizer to: {new_name} (based on first unmatched id)")
        current_normalizer = new_func

        # Remove previous temp file if not the original input
        if temp_in_path != input_fasta and os.path.exists(temp_in_path):
            try:
                os.remove(temp_in_path)
            except Exception:
                pass

        temp_in_path = temp_out_path

    out_handle.close()

    if not keep_temp:
        if os.path.exists(temp_out_path):
            try:
                os.remove(temp_out_path)
            except Exception:
                pass

def move_original_fasta(original_fasta, destination_dir):
    os.makedirs(destination_dir, exist_ok=True)
    destination = os.path.join(destination_dir, os.path.basename(original_fasta))
    try:
        shutil.move(original_fasta, destination)
        print(f"[INFO] Moved original FASTA to: {destination}")
    except Exception as e:
        print(f"[WARNING] Could not move the FASTA file. Error: {e}")

def cleanup_files(paths):
    for p in paths:
        try:
            if p.exists():
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()
        except Exception as e:
            print(f"[WARNING] Cleanup failed for {p}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Cluster a protein FASTA with CD-HIT and annotate each protein header with its cluster ID.")
    parser.add_argument("-i", "--input", help="Input FASTA/FAA file (default: first file in ./database/)", default=None)
    parser.add_argument("-o", "--output", help="Output clustered FASTA file (optional). Defaults to <input>_<identity_percent>percent.fasta", default=None)
    parser.add_argument("-c", "--identity", type=float, default=0.95, help="cd-hit -c sequence identity threshold (0..1). Default: 0.95")
    parser.add_argument("-n", "--word-length", type=int, default=5, help="cd-hit -n word length (default: 5 for proteins >=0.9 identity)")
    parser.add_argument("-T", "--threads", type=int, default=0, help="cd-hit -T number of threads (0 = all CPUs). Default: 0")
    parser.add_argument("-M", "--memory", type=int, default=0, help="cd-hit -M max memory in MB (0 = unlimited). Default: 0")
    parser.add_argument("--cdhit-path", type=str, default=None, help="Path to cd-hit binary (override default)")
    parser.add_argument("--cleanup", type=str2bool, default=True, help="Remove cd-hit outputs after annotating (default: true)")
    parser.add_argument("--keep-cdhit", action="store_true", help="Keep cd-hit produced files even if --cleanup is set (overrides cleanup)")
    parser.add_argument("--move-to", type=str, default=None, help="Optional directory to move the original FASTA after successful processing. By default the input is left in place.")

    args = parser.parse_args()

    input_fasta = args.input or find_default_input()
    input_dir = os.path.dirname(os.path.abspath(input_fasta)) or "."
    prefix = os.path.splitext(os.path.basename(input_fasta))[0]

    # Default output name: input + identity percent (e.g., input_95percent.fasta)
    identity_percent = int(round(args.identity * 100))
    if args.output:
        output_fasta = args.output
    else:
        output_fasta = os.path.join(input_dir, f"{prefix}_{identity_percent}percent.fasta")

    # Existing cluster annotations make CD-HIT unnecessary.
    if is_already_processed(input_fasta):
        print(f"[INFO] File '{input_fasta}' already contains cluster annotations.")
        if os.path.abspath(input_fasta) != os.path.abspath(output_fasta):
            shutil.copy2(input_fasta, output_fasta)
            print(f"[INFO] Copied clustered FASTA to: {output_fasta}")
        return

    # Resolve CD-HIT only when clustering is actually required.
    cdhit_exec = resolve_cdhit_path(args.cdhit_path)

    # Step 1: Sort input FASTA deterministically
    sorted_fasta_path = os.path.join(input_dir, f"{prefix}_sorted.fasta")
    try:
        records = list(SeqIO.parse(input_fasta, "fasta"))
        records.sort(key=lambda r: r.id)  # You can change this to r.description if preferred
        with open(sorted_fasta_path, "w") as sorted_fh:
            SeqIO.write(records, sorted_fh, "fasta")
    except Exception as e:
        print(f"[ERROR] Failed to sort input FASTA: {e}")
        return

    # Run cd-hit
    try:
        out_fa_path, clstr_path = run_cdhit_clustering(
            input_fasta=sorted_fasta_path,
            output_fasta=os.path.join(input_dir, f"{prefix}_cdhit_clustered.fasta"),
            identity=args.identity,
            word_len=args.word_length,
            threads=args.threads,
            memory=args.memory,
            cdhit_path=cdhit_exec
        )
    except subprocess.CalledProcessError:
        # Error already printed in run_cdhit_clustering
        return
    except Exception as e:
        print(f"[ERROR] cd-hit failed: {e}")
        return

    # Parse cluster file and annotate original FASTA
    try:
        cluster_map = parse_clstr_file(clstr_path)
    except Exception as e:
        print(f"[ERROR] Failed to parse cluster file: {e}")
        return

    annotate_fasta_with_clusters(input_fasta, cluster_map, output_fasta)
    print(f"[DONE] Clustered FASTA written to: {output_fasta}")

    # Optionally cleanup cd-hit produced files (unless user requested to keep them)
    if args.cleanup and not args.keep_cdhit:
        tracked = [out_fa_path, clstr_path, Path(sorted_fasta_path)]
        cleanup_files(tracked)
    else:
        if args.keep_cdhit:
            print("[INFO] Keeping cd-hit produced files (user requested --keep-cdhit).")

    # Optionally move the original FASTA after successful processing.
    if args.move_to:
        move_original_fasta(input_fasta, args.move_to)

if __name__ == "__main__":
    main()
