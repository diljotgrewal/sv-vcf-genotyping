import bisect
import collections
import gzip
import os

import click
import pandas as pd
import vcf


class BreakpointDatabase(object):
    def __init__(self, breakpoints, id_col='breakpoint_id'):
        self.positions = collections.defaultdict(list)
        self.break_ids = collections.defaultdict(set)

        for idx, row in breakpoints.iterrows():
            for side in ('1', '2'):
                self.positions[(row['chromosome_' + side], row['strand_' + side])].append(row['position_' + side])
                self.break_ids[(row['chromosome_' + side], row['strand_' + side], row['position_' + side])].add(
                    (row[id_col], side))

        for key in self.positions.keys():
            self.positions[key] = sorted(self.positions[key])

    def query(self, row, extend=0):
        exclusion = row['breakpoint_id']
        # raise Exception(row)
        matched_ids = list()

        for side in ('1', '2'):
            chrom_strand_positions = self.positions[(row['chromosome_' + side], row['strand_' + side])]
            idx = bisect.bisect_left(chrom_strand_positions, row['position_' + side] - extend)
            side_matched_ids = list()

            while idx < len(chrom_strand_positions):
                pos = chrom_strand_positions[idx]
                dist = abs(pos - row['position_' + side])

                if pos >= row['position_' + side] - extend and pos <= row['position_' + side] + extend:
                    for break_id in self.break_ids[(row['chromosome_' + side], row['strand_' + side], pos)]:
                        side_matched_ids.append((break_id, dist))

                if pos > row['position_' + side] + extend:
                    break
                idx += 1
            matched_ids.append(side_matched_ids)

        matched_ids_bypos = list()
        for matched_id_1, dist_1 in matched_ids[0]:
            for matched_id_2, dist_2 in matched_ids[1]:
                if matched_id_1[0] == matched_id_2[0] and matched_id_1[1] != matched_id_2[1]:
                    if not matched_id_1[0] == exclusion:
                        matched_ids_bypos.append((dist_1 + dist_2, matched_id_1[0]))

        ids = set([v[1] for v in matched_ids_bypos] + [exclusion])
        return sorted(ids)


class SvVcfData(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.reader = self._get_reader(filepath)
        self.caller = self._get_caller()

    def _get_caller(self):
        header_infos = self.reader.infos.values()

        if any(["lumpy" in str(v).lower() for v in header_infos]):
            return 'lumpy'
        elif any(["gridss" in str(v).lower() for v in header_infos]):
            return 'gridss'
        elif any(["svaba" in str(v).lower() for v in header_infos]):
            return 'svaba'
        else:
            raise Exception('unknown caller')

    @staticmethod
    def _get_reader(vcf_file):
        return vcf.Reader(filename=vcf_file)

    def _parse_vcf(self):

        for record in self.reader:

            data = {
                'chrom': record.CHROM,
                'pos': record.POS,
                'ref': record.REF,
                'alt': record.ALT,
                'qual': record.QUAL,
                'id': record.ID,
                'filter': record.FILTER,
            }

            # if data['filter'] == []:
            #     data['filter'] = None
            #
            # if data['filter'] is not None:
            #     continue

            info = record.INFO

            for k, v in info.items():
                if isinstance(v, list):
                    v = ';'.join(map(str, v))
                data[k] = v

            for sample in record.samples:
                sample_name = sample.sample
                sample_data = sample.data
                for k, v in sample_data._asdict().items():
                    if isinstance(v, list):
                        v = ';'.join([str(val) for val in v])
                    k = '{}_{}'.format(sample_name, k)
                    data[k] = v

            yield data

    def _group_bnds(self, calls):
        bnds = {}

        for record in calls:
            if self.caller == 'lumpy' and record['SVTYPE'] == 'INV':
                strands = record['STRANDS'].split(';')

                record['STRANDS'] = strands[0]
                yield (record,)
                record['STRANDS'] = strands[1]
                yield (record,)
            elif record['SVTYPE'] == 'BND':
                if 'MATEID' not in record:
                    continue

                if record['MATEID'] in bnds:
                    yield (record, bnds[record['MATEID']])
                    bnds.pop(record['MATEID'])
                else:
                    bnds[record['id']] = record
            else:
                yield record,

        # assert len(bnds) == 0, bnds

    def _get_mates(self, records):
        ends_with_val = {
            'lumpy': '_1',
            'svaba': ':1',
            'gridss': 'h'
        }

        if records[0]['id'].endswith(ends_with_val[self.caller]):
            mate1, mate2 = records
        else:
            mate2, mate1 = records

        return mate1, mate2

    @staticmethod
    def _get_strand_from_alt(alt):
        """
        If the nucleotide comes first, then it is a "+" facing break at
        that site (see element "W" of the VCF4.2 specs, pg 13), otherwise
        it is "-" at that site. Then, if the bracket is ], then the
        partner breakpoint is "+", otherwise it is left facing
        :param alt:
        :type alt:
        :return:
        :rtype:
        """
        return alt[0].orientation

    def _get_strands(self, mate1, mate2):
        if self.caller == 'lumpy':
            strands_1 = mate1['STRANDS'].split(':')[0]
            strands_2 = mate2['STRANDS'].split(':')[0]
            assert strands_1 == strands_2[::-1]
            return strands_1[0], strands_2[0]
        else:
            strand_1 = self._get_strand_from_alt(mate1['alt'])
            strand_2 = self._get_strand_from_alt(mate2['alt'])

        return strand_1, strand_2

    @staticmethod
    def _process_lumpy_unmatched_record(record):
        record = record[0]

        strands = record['STRANDS'].split(':')[0]
        assert len(strands) == 2

        outdata = {
            'chromosome_1': record['chrom'],
            'position_1': record['pos'],
            'chromosome_2': record['chrom'],
            'position_2': record['END'],
            'strand_1': strands[0],
            'strand_2': strands[1],
            'type': record['SVTYPE']
        }

        return outdata

    def _process_bnd_call(self, record):

        assert len(record) == 2

        mate1, mate2 = self._get_mates(record)
        strand_1, strand_2 = self._get_strands(mate1, mate2)

        assert mate1['SVTYPE'] == mate2['SVTYPE']

        outdata = {
            'chromosome_1': mate1['chrom'],
            'position_1': mate1['pos'],
            'strand_1': strand_1,
            'chromosome_2': mate2['chrom'],
            'position_2': mate2['pos'],
            'strand_2': strand_2,
            'type': mate1['SVTYPE']
        }

        return outdata

    def _filter_low_qual_calls(self, calls):

        for call in calls:

            if len(call) == 1 and self.caller == 'lumpy':

                if call[0]['filter'] and 'LOW_QUAL' in call[0]['filter']:
                    continue
            else:
                assert len(call) == 2

                if call[0]['filter'] and 'LOW_QUAL' in call[0]['filter'] and 'LOW_QUAL' in call[1]['filter']:
                    continue

            yield call

    def fetch(self):
        records = self._parse_vcf()
        records = self._group_bnds(records)
        records = self._filter_low_qual_calls(records)

        for record in records:
            if len(record) == 1 and self.caller == 'lumpy':
                yield self._process_lumpy_unmatched_record(record)
            else:
                yield self._process_bnd_call(record)

    def as_data_frame(self):
        data = [record for record in self.fetch()]

        data = pd.DataFrame(data)
        data['caller'] = self.caller

        data['breakpoint_id'] = data.index

        data['breakpoint_id'] = data['breakpoint_id'].astype(str) + '_' + data['caller']

        return data


def read_destruct(destruct_calls):
    df = pd.read_csv(destruct_calls, sep='\t', dtype={'chromosome_1': str, 'chromosome_2': str})

    colnames = [
        'prediction_id', 'chromosome_1', 'position_1', 'strand_1',
        'chromosome_2', 'position_2', 'strand_2', 'type'
    ]
    df = df[colnames]

    df['breakpoint_id'] = df.prediction_id
    del df['prediction_id']
    df['caller'] = 'destruct'

    df['breakpoint_id'] = df['breakpoint_id'].astype(str) + '_' + df['caller']

    return df


def check_common(x, df_db, calls):
    val = df_db.query(x, extend=500)

    val = sorted(val)

    if len(val) == 1:
        return

    if val[0] not in calls:
        calls[val[0]] = set()

    for v in val[1:]:
        calls[val[0]].add(v)


def get_common_calls(df, df_db):
    calls = {}

    for i, row in df.iterrows():
        check_common(row, df_db, calls)

    new_groups = {}
    for i, (key, vals) in enumerate(calls.items()):
        new_groups[key] = i
        for val in vals:
            new_groups[val] = i

    return new_groups


def write_vcf_output(vcf_file, ids, output_file):
    input_open = gzip.open if vcf_file.endswith('gz') else open
    output_open = gzip.open if output_file.endswith('gz') else open

    with output_open(output_file, 'wt') as writer, input_open(vcf_file, 'rt') as reader:
        for line in reader:
            if line.startswith('#'):
                writer.write(line)
                continue

            brk_id = line.strip().split('\t')[2].split('_')[0]

            if brk_id in ids:
                writer.write(line)


def main(
        destruct_csv, lumpy_vcf, gridss_vcf, svaba_vcf, consensus_vcf_lumpy, consensus_vcf_svaba,
        consensus_vcf_gridss
):
    allcalls = [
        read_destruct(destruct_csv),
        SvVcfData(lumpy_vcf).as_data_frame(),
        SvVcfData(gridss_vcf).as_data_frame(),
        SvVcfData(svaba_vcf).as_data_frame(),
    ]

    allcalls = pd.concat(allcalls)

    allcalls = allcalls[allcalls['chromosome_1'] == '1']

    allcalls_db = BreakpointDatabase(allcalls)

    groups = get_common_calls(allcalls, allcalls_db)

    allcalls['grouped_breakpoint_id'] = allcalls['breakpoint_id'].apply(lambda x: groups.get(x, float("nan")))

    allcalls = allcalls[~ pd.isnull(allcalls.grouped_breakpoint_id)]

    allcalls = allcalls.groupby('grouped_breakpoint_id')

    ids = {'svaba': [], 'lumpy': [], 'gridss': []}

    for _, brkgrp in allcalls:

        # filter multiple calls by same tool in the window
        # without confirmation from another tool
        if len(brkgrp.caller.unique()) == 1:
            continue

        brkgrp['caller'] = ','.join(list(brkgrp['caller']))

        tracked = False
        for v in list(brkgrp['breakpoint_id']):
            if 'lumpy' in v:
                ids['lumpy'].append(v.replace('_lumpy', ''))
                tracked = True
            elif 'svaba' in v:
                ids['svaba'].append(v.replace('_svaba', ''))
                tracked = True
            elif 'gridss' in v:
                ids['gridss'].append(v.replace('_gridss', ''))
                tracked = True

        assert tracked is True

    write_vcf_output(lumpy_vcf, ids['lumpy'], consensus_vcf_lumpy)
    write_vcf_output(svaba_vcf, ids['svaba'], consensus_vcf_svaba)
    write_vcf_output(gridss_vcf, ids['gridss'], consensus_vcf_gridss)


@click.group()
def cli():
    pass


@cli.command()
@click.option('--destruct_csv', required=True, help='CSV file path')
@click.option('--lumpy_vcf', required=True, help='vcf')
@click.option('--svaba_vcf', required=True, help='vcf')
@click.option('--gridss_vcf', required=True, help='vcf')
@click.option('--outdir', required=True, help='output dir')
def run(
        destruct_csv,
        lumpy_vcf,
        svaba_vcf,
        gridss_vcf,
        outdir
):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    consensus_lumpy = os.path.join(outdir, 'lumpy.vcf')
    consensus_svaba = os.path.join(outdir, 'svaba.vcf')
    consensus_gridss = os.path.join(outdir, 'gridss.vcf')

    main(
        destruct_csv, lumpy_vcf, svaba_vcf, gridss_vcf,
        consensus_lumpy, consensus_svaba, consensus_gridss
    )


if __name__ == "__main__":
    cli()
