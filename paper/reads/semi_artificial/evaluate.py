# import math
import os
import subprocess

from pathlib import Path


def create_descendants(output_basedir, count=10,
                       distances=[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                                  0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                  0.8, 0.85, 0.9, 0.95, 1]):
    for i in range(0, count):
        coverage = 1/2**i
        for distance in distances:
            out_dir = output_basedir + '/' + \
                str(coverage) + '/' + str(distance) + '/'
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            for j in range(0, 10):
                subprocess.call([
                    '../../../build/generator',
                    '-m', str(distance), '-f', '../assembled-ecoli/536.fasta',
                    # '-i', '0.01',
                    '-o', out_dir + '536_desc_' + str(j) + '.fasta'])


def split_fasta_files(input_basedir, output_basedir):
    for path in [path for path in Path(input_basedir).glob('**/*.fasta')
                 if path.is_file()]:
        print(path)
        with open(path, 'r') as file:
            contents = file.readlines()
        split_idx = 1
        while contents[split_idx][0] != '>':
            split_idx += 1
        components = str(os.path.splitext(path)[0]).split('/')
        components[0] = output_basedir
        path = '/'.join(components)
        if not os.path.exists(path):
            os.makedirs(path)
        with open(path + '/0.fasta', 'w') as file:
            file.write(''.join(contents[0:split_idx]))
        with open(path + '/1.fasta', 'w') as file:
            file.write(''.join(contents[split_idx:]))


def create_reads(inout_basedir):
    for path in [path for path in Path(inout_basedir).glob('**/*.fasta')
                 if path.is_file()]:
        coverage = str(path).split('/')[2]
        subprocess.call([
            '../../../../art/art_illumina', '-ss', 'HS25', '-sam',
            '-i', str(path), '-l',  '150', '-f', str(coverage),
            '-o', os.path.splitext(path)[0] + '_reads'])


def fastq_to_multi_fasta_files(inout_basedir):
    for input_filename in [path for path in Path(inout_basedir).glob('**/*.fq')
                           if path.is_file()]:
        output_filename = os.path.splitext(input_filename)[0]+'.fasta'
        print(input_filename, '->', output_filename)
        with open(input_filename, 'r') as file:
            lines = file.readlines()

        with open(output_filename, 'w') as file:
            for i in range(0, len(lines), 4):
                file.write('>')
                file.write(lines[i][1:])
                file.write(lines[i+1])


def replace_basedir(path, basedir):
    components = path.split('/')
    for i, c in enumerate(basedir.split('/')):
        components[i] = c
    return '/'.join(components)


def estimate(input_basedir, output_basedir, filename_0, filename_1):
    for path in [f for f in Path(input_basedir).glob('**/' + filename_0)
                 if f.is_file()]:
        basedir = os.path.dirname(path)
        out_filepath = replace_basedir(basedir, output_basedir) + '.dmat'
        if not os.path.exists(os.path.dirname(out_filepath)):
            os.makedirs(os.path.dirname(out_filepath))
        subprocess.call([
            '../../../build/slope-spam',
            '--as-reads', basedir + '/' + filename_0,
            basedir + '/' + filename_1, '-o', out_filepath])


def collect(basedir):
    data = {}
    for input_path in [f for f in Path(basedir).glob('**/*.dmat')
                       if f.is_file()]:
        coverage = str(input_path).split('/')[2]
        distance = str(input_path).split('/')[3]
        data.setdefault(coverage, {})
        data[coverage].setdefault(distance, [])
        with open(input_path, 'r') as file:
            content = file.read()
        v = float(content.split()[3])
        data[coverage][distance].append(v)

    for coverage, v in data.items():
        with open(basedir + '/' + coverage + '.csv', 'w') as file:
            for distance, values in v.items():
                file.write(str(distance) + ' ' + str(distance) +
                           ' ' + ' '.join([str(v) for v in values]) + '\n')


def evaluate():
    create_descendants('simple/descendants')
    split_fasta_files('simple/descendants', 'simple/data')
    create_reads('simple/data')
    fastq_to_multi_fasta_files('simple/data')
    estimate('simple/data', 'simple/assembled_assembled', '0.fasta', '1.fasta')
    collect('simple/assembled_assembled')
    estimate('simple/data', 'simple/assembled_reads',
             '0.fasta', '1_reads.fasta')
    collect('simple/assembled_reads')
    estimate('simple/data', 'simple/reads_reads',
             '0_reads.fasta', '1_reads.fasta')
    collect('simple/reads_reads')


def evaluate_small_distances():
    create_descendants('small_distances/descendants',
                       distances=[0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03,
                                  0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065,
                                  0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1])
    split_fasta_files('small_distances/descendants', 'small_distances/data')
    create_reads('small_distances/data')
    fastq_to_multi_fasta_files('small_distances/data')
    estimate('small_distances/data',
             'small_distances/assembled_assembled', '0.fasta', '1.fasta')
    collect('small_distances/assembled_assembled')
    estimate('small_distances/data', 'small_distances/assembled_reads',
             '0.fasta', '1_reads.fasta')
    collect('small_distances/assembled_reads')
    estimate('small_distances/data', 'small_distances/reads_reads',
             '0_reads.fasta', '1_reads.fasta')
    collect('small_distances/reads_reads')


def evaluate_no_indels():
    create_descendants('no_indels/descendants')
    split_fasta_files('no_indels/descendants', 'no_indels/data')
    create_reads('no_indels/data')
    fastq_to_multi_fasta_files('no_indels/data')
    estimate('no_indels/data', 'no_indels/assembled_assembled',
             '0.fasta', '1.fasta')
    collect('no_indels/assembled_assembled')
    estimate('no_indels/data', 'no_indels/assembled_reads',
             '0.fasta', '1_reads.fasta')
    collect('no_indels/assembled_reads')
    estimate('no_indels/data', 'no_indels/reads_reads',
             '0_reads.fasta', '1_reads.fasta')
    collect('no_indels/reads_reads')


evaluate()
