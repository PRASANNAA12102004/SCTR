#!/usr/bin/env python3

import sys
import os

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}
VALID_BASES = {"A", "C", "G", "T"}   # N is NOT allowed


def parse_fasta(fasta_file):
    seqs = {}
    header = None
    seq_chunks = []

    with open(fasta_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq_chunks).upper()
                header = line[1:].replace(" ", "_")
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header:
            seqs[header] = "".join(seq_chunks).upper()

    return seqs


def check_orf(seq):
    failures = []
    seq_len = len(seq)

    # --- Character checks ---
    invalid_positions = []
    n_positions = []

    for i, b in enumerate(seq):
        if b == "N":
            n_positions.append(i + 1)
        elif b not in VALID_BASES:
            invalid_positions.append((i + 1, b))

    if n_positions:
        failures.append(
            "N bases at positions: " +
            ", ".join(map(str, n_positions))
        )

    if invalid_positions:
        failures.append(
            "Invalid characters at positions: " +
            ", ".join([f"{pos}({base})" for pos, base in invalid_positions])
        )

    # --- Length multiple of 3 ---
    if seq_len % 3 != 0:
        failures.append(
            f"Sequence length {seq_len} not divisible by 3"
        )

    # --- Start codon ---
    if seq_len < 3 or seq[:3] != START_CODON:
        failures.append(
            f"Missing or incorrect start codon at positions 1–3 (found: {seq[:3]})"
        )

    # --- Stop codon ---
    if seq_len < 3 or seq[-3:] not in STOP_CODONS:
        failures.append(
            f"Missing or incorrect stop codon at positions {seq_len-2}–{seq_len} "
            f"(found: {seq[-3:]})"
        )

    # --- Premature stop codons ---
    premature_stops = []
    for i in range(3, seq_len - 3, 3):
        codon = seq[i:i + 3]
        if codon in STOP_CODONS:
            premature_stops.append((i + 1, codon))

    if premature_stops:
        failures.append(
            "Premature stop codons at positions: " +
            ", ".join([f"{pos}({codon})" for pos, codon in premature_stops])
        )

    return failures


def write_fasta(handle, header, seq, width=60):
    handle.write(f">{header}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")


def main():
    if len(sys.argv) != 2:
        sys.stderr.write(f"Usage: {sys.argv[0]} <multi_fasta.fa>\n")
        sys.exit(1)

    fasta_file = sys.argv[1]
    prefix = os.path.splitext(fasta_file)[0]

    out_fasta = f"{prefix}.cleaned.fasta"
    out_qc = f"{prefix}.qc.txt"

    sequences = parse_fasta(fasta_file)

    with open(out_fasta, "w") as fasta_out, open(out_qc, "w") as qc_out:
        for seq_id, seq in sequences.items():
            failures = check_orf(seq)

            if failures:
                qc_out.write(f"Sequence: {seq_id}\n")
                qc_out.write(f"Length: {len(seq)}\n")
                for f in failures:
                    qc_out.write(f"  - {f}\n")
                qc_out.write("\n")
                continue

            # Passed QC → remove terminal stop codon
            trimmed_seq = seq[:-3]
            write_fasta(fasta_out, seq_id, trimmed_seq)

    print(f"Cleaned FASTA written to: {out_fasta}")
    print(f"QC report written to:    {out_qc}")


if __name__ == "__main__":
    main()
