from wgs_stats import Gene
import pandas as pd, numpy as np
import argparse, math
from scipy.spatial import ConvexHull, distance
from scipy.spatial.distance import cdist 
from pyproj import Geod
from shapely.geometry import Polygon, polygon

def parse_args():
    parser = argparse.ArgumentParser(description="Record spatial area and locations occupied by variants within a given gene")
    parser.add_argument('--gene', type=str, help='Gene name to analyze', required=True)
    parser.add_argument('--transcript', type=str, help='Transcript ID from gene to analyze (optional)', default=None)
    parser.add_argument('--vcf', type=str, help='Path to the VCF file', required=True)
    parser.add_argument('--filter', type=str, help='Path to the filter VCF file (optional)', default=None)
    parser.add_argument('--annotation', help='Path to genome annotation DB (created by parse_gff.py)', type=str, required=True)
    parser.add_argument('--fasta', type=str, help='Path to the reference genome FASTA file', required=True)
    parser.add_argument('--proteinfasta', type=str, help='Path to the reference protein FASTA file (for checking amino acid sequences)', default=None)
    parser.add_argument('--metadata', type=str, help='Path to metadata file with (longitude, latitude) or (x, y) columns', default=None, required=True)
    parser.add_argument('--longlat', default=True, action='store_false', help='If set, indicates that the coordinates are in a projected coordinate system (x, y) rather than (longitude, latitude)')
    parser.add_argument('--minlocs', type=int, default=3, help='Minimum number of unique locations required to calculate area and max distance (default: 3)')
    parser.add_argument('--savelocs', default=False, action='store_true', help='If set, saves the locations of each variant to a separate file')
    parser.add_argument('--output', type=str, help='Output prefix (appended with genename-genetranscript)', default='./')
    return parser.parse_args()

def polygon_area(coords):
    """
    calculate area over a set of coordinates in square KM
    """
    poly = Polygon(coords)
    poly = polygon.orient(poly)
    geod = Geod(ellps="WGS84")
    poly_area, poly_perimiter = geod.geometry_area_perimeter(poly)
    return poly_area/1000000

def annotate_variant(gene, transcript,
                     position, carrier_idxs, carrier_chrs, 
                     cds_positions, varcds, varcdsidx, varspdxidx, varcodonidx):
    if position < min(cds_positions):
        return 'upstream'
    elif position > max(cds_positions):
        return 'downstream'
    elif position in varcds:
        codonidx = varcodonidx[np.where(varcds == position)[0][0]]
        derivedA = gene.transcripts[transcript]['sequences']['aa_seq_array'][codonidx, carrier_idxs[0]][carrier_chrs[0]]
        ancestrA = gene.transcripts[transcript]['sequences']['wt_aa_seq'][codonidx]
        if derivedA == ancestrA:
            return 'synonymous'
        elif codonidx == 0:
            return 'start_lost'
        elif derivedA == '*' and ancestrA != '*':
            return 'nonsense'
        elif derivedA != '*' and ancestrA == '*':
            return 'stop_lost'
        else:
            return 'nonsynonymous'
    else:
        return 'intronic'

def get_loc_frequencies(gene, transcript, idx, carrier_idxs, carrier_locs, alt):
    locsx = []
    locsy = []
    freqs = []
    nsamples = []
    for loc in carrier_locs:
        residents = np.where(gene.locs == loc)[0]
        genotypes = gene.transcripts[transcript]['sequences']['genotype_array'][idx][residents]
        freq = np.sum(genotypes == alt) / (len(residents) * 2)
        if freq > 0:
            locsx.append(loc[0])
            locsy.append(loc[1])
            freqs.append(freq)
            nsamples.append(len(residents))
    return np.array(locsx), np.array(locsy), freqs, nsamples

def get_variant_stats(gene, transcript, minlocs=3, savelocs=False):
    positions = []
    alt_alleles = []
    annotations = []
    num_carriers = []
    num_locations = []
    frequencies = []
    area_km2 = []
    max_dist_km = []
    if savelocs:
        Lpositions = []
        Lalt_alleles = []
        Lannotations = []
        Lfrequencies = []
        Llocsx = []
        Llocsy = []
        Lnsamples = []
    for idx, site in enumerate(gene.transcripts[transcript]['sequences']['genotype_array']):
        position = gene.transcripts[transcript]['sequences']['positions'][idx]

        # figure out mutation types (synonymous, nonsynonymous, etc)
        # array of genome coding positions
        if gene.transcripts[transcript]['strand'] == '-':
            cds_positions = np.concatenate([np.flip(np.arange(s, e+1)) for s, e in zip(gene.transcripts[transcript]['cds_start'], gene.transcripts[transcript]['cds_end'])])
        else:
            cds_positions = np.concatenate([np.arange(s, e+1) for s, e in zip(gene.transcripts[transcript]['cds_start'], gene.transcripts[transcript]['cds_end'])])
        # variant coding sequence positions, their indices in the coding sequence array, and their indices in the variant array
        varcds, varcdsidx, varspdxidx = np.intersect1d(cds_positions, gene.transcripts[transcript]['sequences']['positions'], return_indices=True)
        varcodonidx = np.floor(varcdsidx / 3).astype(int)
        for alt in np.unique(site):
            if alt < 1:
                continue
            carrier_idxs, carrier_chrs = np.where(site == alt)
            carrier_locs = np.unique(gene.locs[carrier_idxs], axis=0)
            if len(carrier_locs) < minlocs:
                continue
            variant_ann = annotate_variant(gene, transcript,
                                       position, carrier_idxs, carrier_chrs,
                                       cds_positions, varcds, varcdsidx, varspdxidx, varcodonidx)
            variant_freq = np.sum(site.flatten() == alt) / len(site.flatten())
            positions.append(position)
            alt_alleles.append(alt)
            frequencies.append(variant_freq)
            annotations.append(variant_ann)
            num_carriers.append(len(carrier_idxs))
            num_locations.append(len(carrier_locs))
            if savelocs:
                locsx, locsy, freqs, nsamples = get_loc_frequencies(gene, transcript, idx, carrier_idxs, carrier_locs, alt)
                #print(position, alt, carrier_idxs, carrier_locs, len(locsx))
                Lpositions.extend([position]*len(locsx))
                Lalt_alleles.extend([alt]*len(locsx))
                Lannotations.extend([variant_ann]*len(locsx))
                Lfrequencies.extend(freqs)
                Llocsx.extend(locsx)
                Llocsy.extend(locsy)
                Lnsamples.extend(nsamples)
            try:
                hull = ConvexHull(carrier_locs)
                hull_points = carrier_locs[hull.vertices]
                area_km2.append(polygon_area(hull_points))
                dists = cdist(hull_points, hull_points, metric='euclidean')
                max_dist = np.max(dists)
                max_dist_km.append(max_dist * 111) # rough conversion from degrees to KM
            except:
                area_km2.append(math.nan)
                max_dist_km.append(math.nan)
            continue
    stats = {'position': positions,
             'alt_allele': alt_alleles,
             'annotation': annotations,
             'num_carriers': num_carriers,
             'num_locations': num_locations,
             'area_km2': area_km2,
             'max_dist_km': max_dist_km,
             'frequency': frequencies}
    if savelocs:
        locstats = {'position': Lpositions,
                    'alt_allele': Lalt_alleles,
                    'annotation': Lannotations,
                    'frequency': Lfrequencies,
                    'locx': Llocsx,
                    'locy': Llocsy,
                    'nsamples': Lnsamples}
    else:
        locstats = None
    stats = zscore(pd.DataFrame(stats)).to_dict(orient='list')
    return stats, locstats

def zscore(df):
    df = df.sort_values(by=['frequency']).reset_index(drop=True)
    zscore_func = lambda x: (x - x.mean()) / x.std()
    df.insert(2, 'zscore_area', df.groupby(df.index // round(len(df) / 200))['area_km2'].transform(zscore_func)) # ~200 snps per bin
    return df

def main():
    args = parse_args()
    metadata = pd.read_csv(args.metadata, sep='\t')
    gene = Gene(name = args.gene)
    gene.fetch_gene_coordinates(args.annotation, args.fasta)
    gene.fetch_variation(vcf_file = args.vcf,
                         annotation = args.annotation,
                         fasta = args.fasta,
                         filter = args.filter, 
                         samples=metadata,
                         locs=True)
    if args.transcript is None:
        for transcript in gene.transcripts.keys():
            stats, locstats = get_variant_stats(gene, transcript, savelocs=args.savelocs)
            stats['gene'] = [args.gene]*len(stats['position'])
            stats['chromosome'] = [gene.chromosome]*len(stats['position'])
            stats['transcript'] = [transcript]*len(stats['position'])
            pd.DataFrame(stats).to_csv(f"{args.output}{args.gene}_{transcript}_spatial_scan.txt", sep='\t')
            pd.DataFrame(locstats).to_csv(f"{args.output}{args.gene}_{transcript}_spatial_scan_locs.txt", sep='\t') if args.savelocs and locstats is not None else None
    else:
        transcript = args.transcript
        stats, locstats = get_variant_stats(gene, transcript, savelocs=args.savelocs)
        stats['gene'] = [args.gene]*len(stats['position'])
        stats['chromosome'] = [gene.chromosome]*len(stats['position'])
        stats['transcript'] = [transcript]*len(stats['position'])
        pd.DataFrame(stats).to_csv(f"{args.output}{args.gene}_{transcript}_spatial_scan.txt", sep='\t')
        pd.DataFrame(locstats).to_csv(f"{args.output}{args.gene}_{transcript}_spatial_scan_locs.txt", sep='\t') if args.savelocs and locstats is not None else None
    #print(stats)

if __name__ == "__main__":
    main()
    