#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
import datetime
import sys
import edlib

from collections import deque


__all__ = ["deduplicate_bam_exact", "deduplicate_bam_with_errors"]


def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def deduplicate_bam_exact(
    bam_input: str,
    bam_output: str,
    *,
    max_distance: int = 5,          # ends bp window for UMI overlap
    flush_step: int = 1_000_000,    # positional flush interval
    log_every: int = 10_000_000,    # logging interval
    threads: int = 8,               # worker threads
) -> None:
    """
    Parameters
    ----------
    bam_input : str
        Path to the input BAM.
    bam_output : str
        Path to the output BAM.
    max_distance : int, default=5
        Ends bp window for UMI overlap.
    flush_step : int, default=1_000_000
        Positional flush interval.
    log_every : int, default=10_000_000
        Emit a progress log every N reads/records.
    threads : int, default=6
        Worker threads for I/O/compute where applicable.
    """

    # === Open input/output BAMs ===
    bam_in = pysam.AlignmentFile(bam_input, "rb", threads=threads)
    bam_out = pysam.AlignmentFile(bam_output, "wb", header=bam_in.header, threads=threads)
    
    seen = {}  # {(CR, UB): (ref_name, start_pos, end_pos, best_read, count)}
    current_ref = None
    last_flush_pos = 0
    read_count = 0
    output_count = 0
    skipped_count = 0

    print("[Main] Starting processing...")
    
    # Process each contig separately to enable threading
    for ref in bam_in.references:
        print(f"[Main] Processing contig: {ref}")
        for read in bam_in.fetch(ref):
            read_count += 1
            if read.is_unmapped:
                continue
            try:
                cr = read.get_tag("CB")
                ub = read.get_tag("UB")
            except KeyError:
                skipped_count += 1
                continue  # Skip reads without CR or UB
    
            strand = read.is_reverse
            key = (cr, ub, strand)
    
            ref_name = read.reference_name
            start_pos = read.reference_start
            end_pos = read.reference_end
    
            if current_ref != ref_name or start_pos > last_flush_pos + flush_step:
                # Flush out-of-range reads
                keys_to_flush = []
                for k, (rname, spos, epos, best_read, duplicates_count) in seen.items():
                    if rname != ref_name or epos < start_pos:
                        best_read.set_tag(tag="DC", value=duplicates_count, value_type="i")
                        bam_out.write(best_read)
                        output_count += 1
                        keys_to_flush.append(k)
                for k in keys_to_flush:
                    del seen[k]
                current_ref = ref_name
                last_flush_pos = start_pos
    
            if key in seen:
                rname, spos, epos, best_read, duplicates_count = seen[key]
                if rname == ref_name and (
                    (spos - max_distance) <= start_pos <= (spos + max_distance)
                    or (epos - max_distance) <= end_pos <= (epos + max_distance)
                ):
                    if read.query_length > best_read.query_length:
                        seen[key] = (ref_name, start_pos, end_pos, read, duplicates_count + 1)
                    else:
                        seen[key] = (rname, spos, epos, best_read, duplicates_count + 1)
                else:
                    best_read.set_tag(tag="DC", value=duplicates_count, value_type="i")
                    bam_out.write(best_read)
                    output_count += 1
                    seen[key] = (ref_name, start_pos, end_pos, read, 1)
            else:
                seen[key] = (ref_name, start_pos, end_pos, read, 1)
    
            if read_count % log_every == 0:
                print(f"[{now()}] [Main] Processed {read_count:,} reads...")
    
    # Final flush
    for _, (_, _, _, best_read, duplicates_count) in seen.items():
        best_read.set_tag(tag="DC", value=duplicates_count, value_type="i")
        bam_out.write(best_read)
        output_count += 1
    
    bam_in.close()
    bam_out.close()
    
    print(f"[{now()}] [Main] Done.")
    print(f"[Main] Total reads processed: {read_count:,}")
    print(f"[Main] Reads with missing CB/UB tags: {skipped_count:,}")
    print(f"[Main] Final deduplicated reads written: {output_count:,}")


# === Edit distance using edlib ===
def umi_distance(umi1, umi2):
    if len(umi1) != len(umi2):
        return max_edit_distance + 1
    result = edlib.align(umi1, umi2, mode="NW", task="distance", k=max_edit_distance)
    return result["editDistance"]


# === Deduplicate chunk based on fuzzy UMI, duplicate count, and position ===
def deduplicate_chunk(reads_chunk, bam_out):
    deduped_count = 0

    for cb, reads in reads_chunk.items():
        # Sort by descending DC count
        reads = sorted(reads, key=lambda x: -x[0])

        merged = [False] * len(reads)
        result_clusters = []

        for i in range(len(reads)):
            if merged[i]:
                continue
            dc_i, ub_i, start_i, end_i, read_i = reads[i]
            new_cluster = [(dc_i, ub_i, start_i, end_i, read_i)]
            merged[i] = True

            for j in range(i + 1, len(reads)):
                if merged[j]:
                    continue
                dc_j, ub_j, start_j, end_j, read_j = reads[j]

                if (
                    abs(start_i - start_j) <= max_distance or
                    abs(end_i - end_j) <= max_distance
                ) and read_i.is_reverse == read_j.is_reverse and umi_distance(ub_i, ub_j) <= max_edit_distance:
                    new_cluster.append((dc_j, ub_j, start_j, end_j, read_j))
                    merged[j] = True

            result_clusters.append(new_cluster)

        for cluster in result_clusters:
            total_dc = sum(dc for dc, _, _, _, _ in cluster)
            # Use longest read as representative
            best_read = max(cluster, key=lambda x: x[4].query_length)[4]
            # Use UMI from the read with the highest DC
            best_umi = max(cluster, key=lambda x: x[0])[1]
            best_read.set_tag("UB", best_umi)
            best_read.set_tag("DC", total_dc, value_type="i")
            bam_out.write(best_read)
            deduped_count += 1

    return deduped_count


def deduplicate_bam_with_errors(
    bam_input: str,
    bam_output: str,
    *,
    max_distance: int = 5,        # bp for start/end overlap
    max_edit_distance: int = 1,   # max allowed edit distance between UMIs
    log_every: int = 200_000,     # progress log interval
    threads: int = 8,             # threads for I/O
) -> None:
    """
    Parameters
    ----------
    bam_input : str
        Path to the input BAM.
    bam_output : str
        Path to the output BAM.
    max_distance : int, default=5
        bp window for start/end overlap.
    max_edit_distance : int, default=1
        Maximum allowed edit distance between UMIs.
    log_every : int, default=200_000
        Emit a progress log every N records.
    threads : int, default=8
        Threads for I/O/compute where applicable.
    """
	
	# === Main execution ===
    bam_in = pysam.AlignmentFile(bam_input, "rb", threads=threads)
    bam_out = pysam.AlignmentFile(bam_output, "wb", header=bam_in.header, threads=threads)
    
    print(f"[{now()}] [Fuzzy Dedup] Starting second-pass fuzzy UMI deduplication...")
    
    total_input = 0
    total_output = 0
    skipped = 0
    
    for ref in bam_in.references:
        print(f"[{now()}] [Fuzzy Dedup] Processing contig: {ref}")
        current_reads = {}
        max_end_position = 0
        contig_count = 0
    
        for read in bam_in.fetch(ref):
            contig_count += 1
            if contig_count % log_every == 0:
                print(f"[{now()}] [{ref}] Parsed {contig_count:,} reads on contig.")
    
            if read.is_unmapped:
                continue
            try:
                cr = read.get_tag("CB")
                ub = read.get_tag("UB")
                duplicate_count = read.get_tag("DC")
            except KeyError:
                skipped += 1
                continue
    
            start = read.reference_start
            end = read.reference_end
            total_input += 1
    
            if start > max_end_position:
                # Flush current region
                total_output += deduplicate_chunk(current_reads, bam_out)
                current_reads = {}
    
            if cr not in current_reads:
                current_reads[cr] = []
            current_reads[cr].append((duplicate_count, ub, start, end, read))
            max_end_position = max(max_end_position, end)
    
        # Final flush for last chunk in contig
        if current_reads:
            total_output += deduplicate_chunk(current_reads, bam_out)
    
    bam_in.close()
    bam_out.close()
    
    print(f"[{now()}] [Fuzzy Dedup] Done.")
    print(f"[{now()}] [Fuzzy Dedup] Total reads parsed: {total_input:,}")
    print(f"[{now()}] [Fuzzy Dedup] Skipped (missing CB/UB): {skipped:,}")
    print(f"[{now()}] [Fuzzy Dedup] Final deduplicated reads written: {total_output:,}")

