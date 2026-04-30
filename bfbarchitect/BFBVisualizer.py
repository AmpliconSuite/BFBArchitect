import matplotlib.pyplot as plt
import argparse
import pandas as pd
from pathlib import Path
import hashlib

try:
    from bfbarchitect.datatypes import SV, CHR_CENTRO, build_centromere_dict
    from bfbarchitect.utils import create_logger
except:
    from datatypes import SV, CHR_CENTRO, build_centromere_dict
    from utils import create_logger

font = {'family' : 'Arial',
        'size'   : 12}
plt.rc('font', **font)
plt.rcParams['pdf.fonttype'] = 42


def parse_segment_coordinates(file_dir, seg_num):
    segments_coordinates = {}
    with open(file_dir, 'r') as f:
        cnt = 0    
        for line in f:
            if line.startswith('sequence'):
                line = line.strip().split('\t')
                name = chr(ord('A') + cnt)
                cnt += 1
                chrom =  line[1].split(':')[0]
                # if chrom != 'chr13':
                #     continue
                start = int(line[1].split(':')[1][:-1])
                end = int(line[2].split(':')[1][:-1])
                # if not (107462692 <= start and end <= 107617365):
                #     continue
                cn = float(line[3])
                segments_coordinates[name] = {'chrom':chrom, 'start':start, 'end':end,'cn':cn}
                if len(segments_coordinates) == seg_num:
                    break
    return segments_coordinates

def parse_scores(file_dir, multiple = False):
    ans, seg_num = [], 0
    with open(file_dir, 'r') as f: 
        for line in f:
            if len(ans) >= 1 and multiple == False:
                break
            if line.startswith('Path'):
                # create a dictionary to store line information
                info = {}
                for item in line.strip().split(';'):
                    key, value = item.split('=')
                    info[key] = value
                final_score = info['Score']
                structure = ''
                for seg in info['Segments'].split(','):
                    structure += chr(ord('A') + int(seg[:-1]) - 1)
                    seg_num = max(seg_num, abs(int(seg[:-1])))
                scores = {'Structure':structure,'Multiplicity':info.get('Multiplicity', 1),'Final_score':final_score}
                ans.append(scores)
    return ans, seg_num

def parse_foldback_coordinate(file_dir, segments, chrom, start, end):
    foldbacks = []
    # print(chrom, start, end)
    with open(file_dir, 'r') as f:
        for line in f:
            if line.startswith('discordant'):
                line = line.strip().split('\t')
                cn = float(line[2])
                pos1, pos2 = line[1].split('->')
                chr1, chr2 = pos1.split(':')[0], pos2.split(':')[0]
                bp1, bp2 = int(pos1.split(':')[1][:-1]), int(pos2.split(':')[1][:-1])
                str1, str2 = pos1[-1], pos2[-1]
                sv = SV(chr1, bp1, str1, chr2, bp2, str2)
                if sv.type == 'FBI' and sv.chrom1 == chrom and cn > 0:
                    if (start <= sv.bp1 and sv.bp1 <= end) or \
                        (start <= sv.bp2 and sv.bp2 <= end):
                        fb_type = 'left' if sv.strand1 == '-' else 'right'
                        segment_name = None
                        for key in segments.keys():
                            if (fb_type == 'left' and segments[key]['start'] == sv.bp1) or \
                                (fb_type == 'right' and segments[key]['end'] == sv.bp2) :
                                segment_name = key
                                break
                        foldbacks.append({'fb_type':fb_type,'segment_name':segment_name,'start':sv.bp1, 'end':sv.bp2,'chrom':sv.chrom1})
    # print(foldbacks)
    return foldbacks      

def filter_foldbacks(foldbacks,segments_coordinates):
    filtered_fb = []
    for segment in segments_coordinates.keys():
        for fb_type in ['left','right']:
            dist = 9999999999
            ans = ''
            for fb in foldbacks:
                if fb['fb_type'] == fb_type and fb['segment_name'] == segment:
                    avg_point = (fb['start'] + fb['end'])/2
                    if fb_type == 'left':
                        if abs(avg_point - segments_coordinates[segment]['start']) < dist:
                            dist = abs(avg_point - segments_coordinates[segment]['start'])
                            ans = fb
                    if fb_type == 'right':
                        if abs(avg_point - segments_coordinates[segment]['end']) < dist:
                            dist = abs(avg_point - segments_coordinates[segment]['end'])
                            ans = fb
            if ans !='':
                filtered_fb.append(ans)
    return filtered_fb

def detect_start_end(segments):
    chrom = ''
    start = 9999999999
    end = 0
    for k in segments.keys():
        chrom = segments[k]['chrom']
        if segments[k]['start'] < start:
            start = segments[k]['start']
        if segments[k]['end'] > end:
            end = segments[k]['end']
    return chrom, start-300000, end + 300000

def extract_fcna(cnr_fn, chrom, start, end):
    cnr = pd.read_csv(cnr_fn, sep="\t")
    cnr = cnr[(cnr.chromosome == chrom) & (start <= cnr.start) & (cnr.end <= end)]
    x, y = [], []
    x_ranges = []
    for row in cnr.itertuples():
        x.append((row.start+row.end)//2)
        x_ranges.append((row.start, row.end))
        y.append(2**row.log2*2-1)
    return x, x_ranges, y

def plot_segments(ax, segments_coordinates):
    try:
        d = args.deletion.split(',')
        deletion_start, deletion_end = int(d[0]), int(d[1])
    except:
        deletion_start, deletion_end = -1, -1
    for s in segments_coordinates.values():
        if s['start'] <= deletion_start and deletion_end <= s['end']:
            ax.hlines(y = s['cn'], xmin = s['start'], xmax = deletion_start,alpha=0.7,color='black', linewidth=1)
            ax.hlines(y = s['cn'], xmin = deletion_start, xmax = deletion_end,alpha=0.7,color='red', linewidth=1, linestyles='dotted')
            ax.hlines(y = s['cn'], xmin = deletion_end, xmax = s['end'],alpha=0.7,color='black', linewidth=1)
        else:
            ax.hlines(y = s['cn'], xmin = s['start'], xmax = s['end'],alpha=0.7,color='black', linewidth=1)
    for i, s in enumerate(list(segments_coordinates.keys())):
        prev_s = list(segments_coordinates.keys())[i-1]
        if i >0:
            x1 = segments_coordinates[prev_s]['end']
            y2 = segments_coordinates[prev_s]['cn']
            x2 = segments_coordinates[s]['start']
            y1 = segments_coordinates[s]['cn']
            ax.vlines(ymin =min(y1,y2),ymax = max(y1,y2), x = x1,alpha=0.7,color='black', linewidth=1)

def plot_genes(ax, gene_annotation, chrom, start, end, max_y):
    # read genes from gtf file
    oncogenes = dict()
    fp = open(gene_annotation, 'r')
    for line in fp:
        s = line.strip().split('\t')
        if "chr" not in s[0]:
            s[0] = "chr" + s[0]
        if s[0] != chrom:
            continue
        gene_start, gene_end = int(s[3]), int(s[4])
        if gene_end < start or gene_start > end:
            continue
        gene_name = ""
        for token in s[-1].split(';'):
            if "Name" in token:
                gene_name = token[5:]
                break
            if "gene_name" in token:
                gene_name = token.strip()[11:-1]
                break
        if gene_name not in oncogenes:
            oncogenes[gene_name] = [gene_start, gene_end, s[6]]
        else:
            oncogenes[gene_name][0] = min(gene_start, oncogenes[gene_name][0])
            oncogenes[gene_name][1] = max(gene_end, oncogenes[gene_name][1])
    fp.close()
    # plot genes
    gene_colors = dict()
    for gene in oncogenes.keys():
        color = str(hashlib.sha1(gene.encode('utf-8')).hexdigest())[-6:]
        gene_colors[gene] = f'#{color}'
    for gene_name, (gene_start, gene_end, _) in oncogenes.items():
        ax.hlines(y = max_y*1.15, xmin = gene_start, xmax = gene_end, color=gene_colors[gene_name], linewidth=2)
        ax.annotate(gene_name, xy=((gene_start + gene_end)/2, max_y*1.16), ha='center', fontsize=6, color=gene_colors[gene_name])

def plot_foldbacks(ax, foldbacks_coordinate,max_y,max_x,start_x):
    prop = dict(arrowstyle="-|>,head_width=0.1,head_length=0.17",
            shrinkA=0,shrinkB=0,color = 'darkblue',alpha = 0.7, linewidth = 0.5)
    arrow_length = 0.03 * max_x
    for foldback in foldbacks_coordinate:
        ax.plot([foldback['start'], foldback['end']], [max_y*1.08,max_y*1.13],color='darkblue',linestyle='--', alpha = 0.6, linewidth=0.5)
        if foldback['fb_type'] == 'left':
            ax.annotate("", xy=(foldback['start'],max_y*1.08), xytext=(foldback['start']+arrow_length,max_y*1.08), arrowprops=prop)
            ax.annotate("", xy=(foldback['end']+arrow_length,max_y*1.13), xytext=(foldback['end'],max_y*1.13), arrowprops=prop)
        if foldback['fb_type'] == 'right':
            ax.annotate("", xy=(foldback['start'],max_y*1.08), xytext=(foldback['start']-arrow_length,max_y*1.08), arrowprops=prop)
            ax.annotate("", xy=(foldback['end']-arrow_length,max_y*1.13), xytext=(foldback['end'],max_y*1.13), arrowprops=prop)

def plot_rectangle_plot(ax, segments_coordinates,max_cn):
    for s in segments_coordinates:
        ax.add_patch(plt.Rectangle((segments_coordinates[s]['start'], 0), segments_coordinates[s]['end'] - segments_coordinates[s]['start'],
                                      - 0.07 * max(1.2 * max_cn, max_cn + 3), edgecolor='r', facecolor='none',
                                      clip_on=False, alpha=1))
        x = 1.1 * ((segments_coordinates[s]['start'] + segments_coordinates[s]['end']) / 2)
        ax.annotate(s, xy=(x - 0.15 * (segments_coordinates[s]['end'] - segments_coordinates[s]['start']),
                                    - 0.07 * max(1.2 * max_cn, max_cn + 3)), weight="bold")
def plot_segments_border(ax, segments_coordinates,max_cn,max_y):
    ax.vlines(ymin = 0 , ymax = max_cn, x = list(segments_coordinates.values())[0]['start'], linewidth = 1, alpha = 0.8 ,linestyle='--', color = 'gray')
    for s in segments_coordinates:
        ax.vlines(ymin = 0 , ymax = max_cn, x = segments_coordinates[s]['end'], linewidth = 1, alpha = 0.8 ,linestyle='--', color = 'gray')
        x = (segments_coordinates[s]['start'] + segments_coordinates[s]['end']) / 2
        ax.annotate(s, xy=(x ,0.02 * max(1.2 * max_y, max_y + 3)), weight="bold", ha = 'center')

def find_in_foldback(segment, direction, foldbacks):
    for f in foldbacks:
        if f['segment_name'] == segment and f['fb_type'] == direction:
            return True
    return False

def plot_structure(ax, score, segments_coordinates, arm,max_y,max_x, foldbacks):
    try:
        d = args.deletion.split(',')
        deletion_start, deletion_end = int(d[0]), int(d[1])
    except:
        deletion_start, deletion_end = -1, -1
    structure = score['Structure']
    arrow_length = 0.03 * max_x
    prop = dict(arrowstyle="-|>,head_width=0.05,head_length=0.15",
            shrinkA=0,shrinkB=0,color = 'darkblue',alpha = 0.6, linewidth = 0.5)
    rectangle_width = 0.008
    y_increase = 0.028
    color = '#0072b2'
    edgecolor = 'black'
    direction = 'right'
    prev = ''
    y = max_y* 1.2
    if arm == 'p':
        direction = 'left'
    for i, s in enumerate(structure):
        if s != prev:
            start, end = segments_coordinates[s]['start'], segments_coordinates[s]['end']
            if start <= deletion_start and deletion_end <= end:
                ax.add_patch(
                    plt.Rectangle((segments_coordinates[s]['start'], y), deletion_start - segments_coordinates[s]['start'],
                                rectangle_width* max_y, edgecolor=edgecolor, facecolor=color,
                                clip_on=False, alpha=0.6, linewidth=0.3))
                ax.add_patch(
                    plt.Rectangle((deletion_start, y), deletion_end - deletion_start,
                                rectangle_width* max_y, edgecolor='red', facecolor='red',
                                clip_on=False, alpha=0.6, linewidth=0.6, fill=False, linestyle='dotted'))
                ax.add_patch(
                    plt.Rectangle((deletion_end, y), segments_coordinates[s]['end'] - deletion_end,
                                rectangle_width* max_y, edgecolor=edgecolor, facecolor=color,
                                clip_on=False, alpha=0.6, linewidth=0.3))
            else:
                ax.add_patch(
                    plt.Rectangle((segments_coordinates[s]['start'], y), segments_coordinates[s]['end'] - segments_coordinates[s]['start'],
                                rectangle_width* max_y, edgecolor=edgecolor, facecolor=color,
                                clip_on=False, alpha=0.6, linewidth=0.3))
            prev = s
        else:
            y = y + y_increase * max_y
            ax.add_patch(
                plt.Rectangle((segments_coordinates[s]['start'], y), segments_coordinates[s]['end'] - segments_coordinates[s]['start'],
                              rectangle_width* max_y, edgecolor=edgecolor, facecolor=color,
                              clip_on=False, alpha=0.6, linewidth=0.3))   
            if direction == 'right':
                linestyle = '--'
                color2 = 'red'
                prop = dict(arrowstyle="-|>,head_width=0.05,head_length=0.15",
                    shrinkA=0,shrinkB=0,color = 'red',alpha = 0.6, linewidth = 0.5)
                if find_in_foldback(s,direction,foldbacks):
                    linestyle = '-'
                    color2 = 'darkblue'
                    prop = dict(arrowstyle="-|>,head_width=0.05,head_length=0.15",
                        shrinkA=0,shrinkB=0,color = 'darkblue',alpha = 0.6, linewidth = 0.5)
                direction = 'left'
                x_point = [segments_coordinates[s]['end'],segments_coordinates[s]['end'] + arrow_length,segments_coordinates[s]['end'] + arrow_length]
                y_point = [y - y_increase * max_y + (rectangle_width/2) * max_y, y - y_increase * max_y + (rectangle_width/2) * max_y, y + (rectangle_width/2) * max_y]
                ax.plot(x_point, y_point, color=color2, alpha=0.6, linewidth=0.5, linestyle = linestyle)
                ax.annotate("", xy=(segments_coordinates[s]['end'],y + (rectangle_width/2) * max_y), xytext=(segments_coordinates[s]['end'] + arrow_length,y + (rectangle_width/2) * max_y), arrowprops=prop)
            else:
                linestyle = '--'
                color2 = 'red'
                prop = dict(arrowstyle="-|>,head_width=0.05,head_length=0.15",
            shrinkA=0,shrinkB=0,color = 'red',alpha = 0.6, linewidth = 0.5)
                if find_in_foldback(s,direction,foldbacks):
                    linestyle = '-'
                    color2 = 'darkblue'
                    prop = dict(arrowstyle="-|>,head_width=0.05,head_length=0.15",
            shrinkA=0,shrinkB=0,color = 'darkblue',alpha = 0.6, linewidth = 0.5)
                direction = 'right'
                x_point = [segments_coordinates[s]['start'],segments_coordinates[s]['start'] - arrow_length,segments_coordinates[s]['start'] - arrow_length]
                y_point = [y - y_increase * max_y + (rectangle_width/2) * max_y, y - y_increase * max_y + (rectangle_width/2) * max_y, y + (rectangle_width/2) * max_y]
                ax.plot(x_point, y_point, color=color2, alpha=0.6, linewidth=0.5, linestyle = linestyle)
                ax.annotate("", xy=(segments_coordinates[s]['start'],y + (rectangle_width/2) * max_y), xytext=(segments_coordinates[s]['start'] - arrow_length,y + (rectangle_width/2) * max_y), arrowprops=prop)
    ax.annotate('x'+str(score['Multiplicity']), xy = (ax.get_xlim()[1]-0.08*max_x, 0.90*ax.get_ylim()[1]),weight = 'bold',ha = 'center')

def visualize_BFB(cycle_file, graph_file, cnr_file, output_prefix, gene_annotation=None, deletion=None, pdf=False, multiple=False, centromere_dict=None):
    logger = create_logger('BFBVisualizer', f'{output_prefix}_visualization.log')
    logger.info(f'Command: python {Path(__file__).resolve()} --graph {graph_file} --cycle {cycle_file} --cnr {cnr_file} --output_prefix {output_prefix}' + (f' --deletion {deletion}' if deletion else '')
            + (f' --gene {gene_annotation}' if gene_annotation else '') + (f' --pdf' if pdf else ''))
    all_scores, seg_num = parse_scores(cycle_file, multiple=multiple)
    segments_coordinates = parse_segment_coordinates(graph_file, seg_num)
    reconstructed_structure = ''
    chrom , start , end = detect_start_end(segments_coordinates)
    foldbacks_coordinate = parse_foldback_coordinate(graph_file, segments_coordinates, chrom, start, end)
    if cnr_file == None:
        x = []
        y = []
        x_ranges = []
        for s in segments_coordinates.keys():
            x.append((segments_coordinates[s]['start'] + segments_coordinates[s]['end'])//2)
            x_ranges.append((segments_coordinates[s]['start'], segments_coordinates[s]['end']))
            y.append(segments_coordinates[s]['cn'])
    else:
        x, x_ranges, y = extract_fcna(cnr_file, chrom, start, end)
    if centromere_dict is None:
        centromere_dict = CHR_CENTRO
    arm = ''
    if max(x) < centromere_dict.get(chrom, CHR_CENTRO.get(chrom, 0)):
        arm = 'p'
    else:
        arm = 'q'
    for index_1, scores in enumerate(all_scores):
        plt.clf()
        fig, ax = plt.subplots()
        fig.set_size_inches(4, 3.5)
        plt.scatter(x, y, c ="#0072b2", s= 0.1,alpha = 0.5)
        # plt.stackplot(x, y, color='#d55e00', alpha=0.3)
        for i in range(len(x_ranges)):
            x_start, x_end = x_ranges[i]
            height = y[i]
            width = x_end - x_start
            rect = plt.Rectangle((x_start, 0), width, height, edgecolor=None, facecolor='#d55e00', alpha=0.3)
            ax.add_patch(rect)
            
        plot_segments(ax, segments_coordinates)
        if gene_annotation is not None:
            plot_genes(ax, gene_annotation, chrom, start, end, max(y))
        plot_foldbacks(ax, foldbacks_coordinate,max(y),max(x)-min(x), min(x))
        plot_structure(ax, scores, segments_coordinates,arm,max(y),max(x)-min(x),foldbacks_coordinate)
        plt.xlabel('Position', fontsize=16)
        plt.ylabel('Copy number', fontsize=16)
        ylim = ax.get_ylim()
        plt.ylim(0, ylim[1] * 1.05)
        plot_segments_border(ax, segments_coordinates,ax.get_ylim()[1], max(y))
        ytcik = ax.get_yticks()
        new_ytick = []
        for i in ytcik:
            if i < max(y)*1.1:
                new_ytick.append(i)
        plt.yticks(new_ytick)
        ax.add_patch(plt.Rectangle((ax.get_xlim()[0], max(y)*1.04), ax.get_xlim()[1]-ax.get_xlim()[0],
                                        0.13*max(y), edgecolor='none', facecolor='y',
                                        alpha=0.1))
        plt.subplots_adjust(bottom=0.15, left = 0.175)
        plt.legend(handles=[], loc ="upper left",title="Score = {score}".format(score = round(float(scores['Final_score']), 2)),prop={'size': 4},title_fontsize=6)
        # plt.annotate("Score = {score}".format(score = str(scores['Final_score'])[:4]),xy = (0.1,0.1),xycoords='axes fraction',fontsize=6)
        cell_line = output_prefix.split('/')[-1].split('_')[0]
        plt.title(cell_line+' '+str(chrom)+arm, fontsize=16)
        plt.savefig(output_prefix+'_'+str(index_1+1)+('.pdf' if pdf else '.png'), dpi = 300)
        plt.close()
        print('Saved figure to '+output_prefix+'_'+str(index_1+1)+('.pdf' if pdf else '.png'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-graph", "--graph", help="graph.txt dir", required=True)
    parser.add_argument("-cycle", "--cycle", help="cycles.txt dir", required=True)
    parser.add_argument("-cnr", "--cnr", help="cnr dir", required=True)
    parser.add_argument("-o", "--output_prefix", help="Output_dir", required=True)
    parser.add_argument("-d", "--deletion", help="Deletion, e.g., \"1,10000\"")
    parser.add_argument("-pdf", "--pdf", action='store_true', help="Output pdf format")
    parser.add_argument("-g", "--gene", help="Gene annotation", default=None)
    parser.add_argument("-m", "--multiple", action='store_true', help="Visualize all structures")
    parser.add_argument("--centromere", help="Path to a centromere BED file (hg38 defaults used if not provided).", default=None)
    args = parser.parse_args()
    visualize_BFB(args.cycle, args.graph, args.cnr, args.output_prefix, gene_annotation=args.gene, deletion=args.deletion, pdf=args.pdf, multiple=args.multiple, centromere_dict=build_centromere_dict(args.centromere))