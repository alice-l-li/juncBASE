#! /usr/bin/env python3

import sys
import argparse
import math
from statistics import median, stdev

REGIONSIZE = 10 ** 6 * 2


def parse_args():
    parser = argparse.ArgumentParser(
        description='filters and converts the format of juncBase output to make it usable for downstream analyses'
    )
    parser.add_argument(
        '-j', '--juncbase', required=True,
        help='juncbase table with full event info and inclusion;exclusion counts'
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help='output prefix'
    )
    parser.add_argument(
        '--gtf',
        help='[optional] annotated gtf for improving gene definitions'
    )

    args = parser.parse_args()
    return args


def getregion(val):
    return math.floor(val / REGIONSIZE) * REGIONSIZE


def parse_gtf_line(line):
    line = line.split('\t')
    chrom = line[0]
    start, end = int(line[3]), int(line[4])
    gene = line[8].split('gene_name "')[1].split('"')[0]
    strand = line[6]
    return chrom, strand, start, end, gene


def consolidate_junctions(genetotranscripttoexons, genetochromstrand, chromstrandtoregiontojuncs):
    for gene in genetotranscripttoexons:
        allgenejuncs = set()
        for transcript in genetotranscripttoexons[gene]:
            exons = sorted(genetotranscripttoexons[gene][transcript])
            junctions = [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)]
            allgenejuncs.update(junctions)
        chrom, strand = genetochromstrand[gene]
        for j in allgenejuncs:
            chromstrandtoregiontojuncs[(chrom, strand)][getregion(j[0])][j] = gene
    return chromstrandtoregiontojuncs


def parse_gtf_to_juncs_regions(gtffile):
    chromstrandtoregiontojuncs, chromstrandtoregiontogenes = {}, {}
    genetotranscripttoexons = {}
    genetochromstrand = {}

    for line in open(gtffile):
        if line and line[0] != '#':
            etype = line.split('\t', 3)[2]
            if etype == 'gene' or etype == 'exon':
                chrom, strand, start, end, gene = parse_gtf_line(line)
                if etype == 'gene':
                    if (chrom, strand) not in chromstrandtoregiontogenes:
                        chromstrandtoregiontogenes[(chrom, strand)] = {}
                        chromstrandtoregiontojuncs[(chrom, strand)] = {}
                    for region in range(getregion(start), getregion(end) + 1, REGIONSIZE):
                        if region not in chromstrandtoregiontogenes[(chrom, strand)]:
                            chromstrandtoregiontogenes[(chrom, strand)][region] = set()
                            chromstrandtoregiontojuncs[(chrom, strand)][region] = {}
                        chromstrandtoregiontogenes[(chrom, strand)][region].add((start, end, gene))
                    if gene not in genetotranscripttoexons:
                        genetotranscripttoexons[gene] = {}
                        genetochromstrand[gene] = (chrom, strand)
                elif etype == 'exon':
                    transcript = line.split('transcript_id "')[1].split('"')[0]
                    if transcript not in genetotranscripttoexons[gene]:
                        genetotranscripttoexons[gene][transcript] = []
                    genetotranscripttoexons[gene][transcript].append((start, end))

    chromstrandtoregiontojuncs = consolidate_junctions(
        genetotranscripttoexons, genetochromstrand, chromstrandtoregiontojuncs
    )
    return chromstrandtoregiontogenes, chromstrandtoregiontojuncs


def get_correct_posinfo(event_type, exclusion_junctions, inclusion_junctions, inclusion_exons):
    posinfo = None
    if 'cassette' in event_type:
        posinfo = inclusion_exons
    elif 'intron' in event_type:
        posinfo = exclusion_junctions
    elif 'mutually_exclusive' in event_type:
        posinfo = inclusion_junctions
    elif 'alternative' in event_type:
        posinfo = inclusion_junctions
    elif 'exon' in event_type:
        posinfo = inclusion_exons
    return posinfo


def extract_junc_pos(juncstring):
    alljuncs = [x.split(',') for x in juncstring.split(';') if x != '']
    alljuncs = [x for xs in alljuncs for x in xs]
    alljuncs = [x.split(':')[1].split('-') for x in alljuncs]
    alljuncs = [(int(x[0]), int(x[1])) for x in alljuncs]
    return alljuncs


def get_gene_from_juncs(junctocheck, chrom, strand, chromstrandtoregiontojuncs):
    alljuncs = extract_junc_pos(junctocheck)
    allpos = [x for xs in alljuncs for x in xs]
    annotjunctionstocheck = {}

    if (chrom, strand) in chromstrandtoregiontojuncs:
        for roundedpos in range(getregion(min(allpos)), getregion(max(allpos)) + 1, REGIONSIZE):
            if roundedpos in chromstrandtoregiontojuncs[(chrom, strand)]:
                annotjunctionstocheck.update(chromstrandtoregiontojuncs[(chrom, strand)][roundedpos])

    for j in alljuncs:
        if j in annotjunctionstocheck:
            return annotjunctionstocheck[j]

    return None


def get_overlap_genes(mypos, annotgenestocheck):
    mygenes = []
    for gstart, gend, gname in annotgenestocheck:
        overlap = min((gend, mypos[1])) - max((gstart, mypos[0]))
        if overlap > 0 and overlap > (mypos[1] - mypos[0]) / 2:
            mygenes.append((overlap, gname))
    return mygenes


def expand_strand(strand):
    if strand == '.':
        return ['+', '-']
    else:
        return [strand]


def get_gene_from_region(regiontocheck, chrom, strand, chromstrandtoregiontogenes):
    allregions = extract_junc_pos(regiontocheck)
    mypos = [x for xs in allregions for x in xs]
    mypos = (min(mypos), max(mypos))

    annotgenestocheck = set()
    for strand_i in expand_strand(strand):
        if (chrom, strand_i) not in chromstrandtoregiontogenes:
            continue
        for roundedpos in range(getregion(mypos[0]), getregion(mypos[1]) + 1, REGIONSIZE):
            if roundedpos in chromstrandtoregiontogenes[(chrom, strand_i)]:
                annotgenestocheck.update(chromstrandtoregiontogenes[(chrom, strand_i)][roundedpos])

    mygenes = get_overlap_genes(mypos, annotgenestocheck)
    if len(mygenes) > 0:
        mygenes.sort(reverse=True)
        return mygenes[0][1]
    else:
        return None


def identify_gene_for_event(chromstrandtoregiontojuncs, chromstrandtoregiontogenes, junctocheck, regiontocheck, chrom, strand):
    newgene = None
    if junctocheck != '' and strand != '.':
        newgene = get_gene_from_juncs(junctocheck, chrom, strand, chromstrandtoregiontojuncs)

    if not newgene:
        newgene = get_gene_from_region(regiontocheck, chrom, strand, chromstrandtoregiontogenes)

    if newgene:
        return newgene
    else:
        return 'nogene'


def extract_keys_from_juncbase(jbfile, chromstrandtoregiontojuncs, chromstrandtoregiontogenes):
    header = None
    keytoevents = {}

    for raw_line in open(jbfile):
        line = raw_line.rstrip('\n').split('\t')

        if not line:
            continue

        if line[0][0] != '#':
            event_type = line[1]
            if 'jcn_only' not in event_type and not (event_type == 'intron' and line[0] == 'N'):
                exclusion_junctions, inclusion_junctions, exclusion_exons, inclusion_exons = line[5:9]
                posinfo = get_correct_posinfo(event_type, exclusion_junctions, inclusion_junctions, inclusion_exons)
                if posinfo:
                    chrom, strand = line[3], line[4]
                    if chromstrandtoregiontogenes is not None:
                        newgene = identify_gene_for_event(
                            chromstrandtoregiontojuncs,
                            chromstrandtoregiontogenes,
                            inclusion_junctions,
                            posinfo,
                            chrom,
                            strand
                        )
                        line[2] = newgene

                    info = ';'.join(line[:3] + [strand])
                    key = (info, posinfo)
                    if key not in keytoevents:
                        keytoevents[key] = []
                    keytoevents[key].append(line)
        else:
            header = line

    return header, keytoevents


def safe_int(x):
    x = x.strip()
    if x == "":
        return 0
    return int(x)


def calculate_inc_exc_vals(vals_str):
    inclusion_vals = []
    exclusion_vals = []

    for x in vals_str:
        x = x.strip()

        if x == "":
            exc = 0
            inc = 0
        else:
            parts = x.split(';')
            exc = safe_int(parts[0]) if len(parts) > 0 else 0
            inc = safe_int(parts[1]) if len(parts) > 1 else 0

        inclusion_vals.append(inc)
        exclusion_vals.append(exc)

    tot_counts = [inclusion_vals[i] + exclusion_vals[i] for i in range(len(inclusion_vals))]
    PSI = [inclusion_vals[i] / tot_counts[i] if tot_counts[i] > 0 else 'NA' for i in range(len(inclusion_vals))]
    return inclusion_vals, exclusion_vals, tot_counts, PSI


def safe_stdev(values):
    vals = [x for x in values if x != 'NA']
    if len(vals) < 2:
        return 0
    return stdev(vals)


def basic_dedup_build_groups(keytoevents):
    einfotoepos = {}
    for key in keytoevents:
        mydata = keytoevents[key]

        for i in range(len(mydata)):
            inclusion_vals, exclusion_vals, tot_counts, PSI = calculate_inc_exc_vals(mydata[i][11:])
            mydata[i].extend([inclusion_vals, exclusion_vals, tot_counts, PSI])

        mydata.sort(
            key=lambda x: median(x[-2]) * safe_stdev(x[-1]),
            reverse=True
        )
        mydata = mydata[0]

        einfo, eposinfo = key
        allpos = extract_junc_pos(eposinfo)
        allpos = [x for xs in allpos for x in xs]

        if len([x for x in mydata[-2] if x != 0]) >= 3:
            if einfo not in einfotoepos:
                einfotoepos[einfo] = []
            einfotoepos[einfo].append(((min(allpos), max(allpos)), eposinfo, mydata))

    return einfotoepos


def generate_overlap_groups(eventposlist):
    newgroups, group = [], []
    lastend = 0
    eventposlist.sort()

    for eventposinfo in eventposlist:
        start, end = eventposinfo[0]
        if start > lastend:
            if len(group) > 0:
                newgroups.append(group)
            group = []

        group.append(eventposinfo)
        if end > lastend:
            lastend = end

    if len(group) > 0:
        newgroups.append(group)

    return newgroups


def collapse_groups_write_out(outprefix, header, einfotoepos):
    with open(outprefix + '.drimformat.tsv', 'w') as out, \
         open(outprefix + '.jb.dedup.psi.tsv', 'w') as outpsi, \
         open(outprefix + '.dedup.bed', 'w') as outbed:

        out.write('\t'.join(['', 'gene_id', 'feature_id'] + header[11:]) + '\n')
        outpsi.write('\t'.join(['#drim_id'] + header) + '\n')

        c = 0
        for einfo in einfotoepos:
            newgroups = generate_overlap_groups(einfotoepos[einfo])

            for group in newgroups:
                group.sort(
                    key=lambda x: median(x[-1][-2]) * safe_stdev(x[-1][-1]),
                    reverse=True
                )
                myevent = group[0]
                gene_id = einfo + ';' + myevent[1]
                inc, exc, tot, psi = myevent[-1][-4:]

                outbed.write(
                    '\t'.join([
                        myevent[1].split(':')[0],
                        str(myevent[0][0]),
                        str(myevent[0][1]),
                        gene_id
                    ]) + '\n'
                )

                feature_id = 'inclusion_' + gene_id
                out.write('\t'.join([str(c), gene_id, feature_id] + [str(x) for x in inc]) + '\n')
                c += 1

                feature_id = 'exclusion_' + gene_id
                out.write('\t'.join([str(c), gene_id, feature_id] + [str(x) for x in exc]) + '\n')
                c += 1

                outpsi.write(
                    '\t'.join(
                        [gene_id] + myevent[-1][:11] +
                        [str(round(100 * x, 2)) if x != 'NA' else x for x in psi]
                    ) + '\n'
                )


if __name__ == "__main__":
    args = parse_args()

    chromstrandtoregiontojuncs, chromstrandtoregiontogenes = None, None
    if args.gtf:
        sys.stderr.write('processing gtf\n')
        chromstrandtoregiontogenes, chromstrandtoregiontojuncs = parse_gtf_to_juncs_regions(args.gtf)

    sys.stderr.write('reading in juncbase table\n')
    header, keytoevents = extract_keys_from_juncbase(
        args.juncbase,
        chromstrandtoregiontojuncs,
        chromstrandtoregiontogenes
    )

    sys.stderr.write('building groups and doing basic dedup\n')
    einfotoepos = basic_dedup_build_groups(keytoevents)

    sys.stderr.write('collapsing groups and writing output\n')
    collapse_groups_write_out(args.output, header, einfotoepos)