import os
import subprocess
import yaml


def split(software):
    software = software.strip()
    software = software.replace(' ', '').replace(',', '')
    for i, c in enumerate(software):
        if c in ('<', '>', '=', '!', '~'):
            return software[:i], software[i:]
    raise ValueError('No version condition found in {}'.format(software))


def split_conditions(conditions):
    """
    example:
    '>=3.4.5!=3.5.0<4.0' -> [('>=', '3.4.5'), ('!=', '3.5.0'), ('<', '4.0')]
    """
    out = []
    i = 0
    while i < len(conditions):
        for operator in ('>=', '>', '<=', '<', '!=', '='):
            if conditions.startswith(operator, i):
                i += len(operator)
                j = i
                while j < len(conditions) and not conditions[j] in ('>=', '>', '<=', '<', '!=', '=', '!'):
                    j += 1
                out.append((operator, conditions[i:j]))
                i = j
                break
        else:
            raise ValueError('Invalid operator in conditions')
    return out


def is_number_or_dot(s):
    for c in s:
        if c not in ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.'):
            return False
    return True


def compare_version_numbers(a, b):
    """
    Example:
    '1.2.3', '1.2.3' -> 0
    '1.2.3', '1.2.4' -> -1
    '1.2.4', '1.2.3' -> 1
    """

    # compare numbers
    a = a.split('.')
    b = b.split('.')
    # make sure they are the same length
    while len(a) < len(b):
        a.append('0')
    while len(b) < len(a):
        b.append('0')
    a, b = ''.join(a), ''.join(b)
    if a == b:
        return 0
    return (a > b) - (a < b)


def validate_condition(installed_version, required_version, condition):
    if condition == '=':
        return installed_version == required_version
    if condition == '>':
        return compare_version_numbers(installed_version, required_version) > 0
    if condition == '<':
        return compare_version_numbers(installed_version, required_version) < 0
    if condition == '>=':
        return compare_version_numbers(installed_version, required_version) >= 0
    if condition == '<=':
        return compare_version_numbers(installed_version, required_version) <= 0
    if condition == '!=':
        return installed_version != required_version
    raise Exception('condition {} not recognized'.format(condition))


def validate_all_conditions(installed_version, depends):
    """Example:
    ('fastqc', '0.11.8', [('>=', '0.11.8'), ('!=', '0.11.9')]) -> True
    ('fastqc', '0.11.9', [('>=', '0.11.8'), ('!=', '0.11.9')]) -> False
    """
    for condition, required_version in depends:
        if not validate_condition(installed_version, required_version, condition):
            return False
    return True


dependencies_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'env.yml')

dependencies_yml = yaml.load(open(dependencies_file), Loader=yaml.FullLoader)['dependencies']
del dependencies_yml[-1]

dependencies_txt = {split(v)[0]: split(v)[1] for v in dependencies_yml}
dependencies = {k: split_conditions(split(v)[1]) for k, v in dependencies_txt.items()}
requirements = open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'requirements.txt')).readlines()
dependencies_txt = {**dependencies_txt, **{split(v)[0]: split(v)[1] for v in requirements}}
use_requirements = ['cutadapt']
requirements = {split(v)[0]: split_conditions(split(v)[1]) for v in requirements if split(v)[0] in use_requirements}
dependencies = {**dependencies, **requirements}

dependencies_installation_page = 'https://rouskinlab.github.io/dreem/dreem/dependencies.html'


def run_cmd(cmd):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output.decode('utf-8').strip()


def check_fastqc_exists():
    version = dependencies['fastqc']
    version_txt = dependencies_txt['fastqc']
    fastqc_installed_version = run_cmd(['fastqc', '--version'])
    assert fastqc_installed_version.startswith('FastQC'), 'fastqc is not installed'
    fastqc_installed_version = fastqc_installed_version.split('\n')[0].split(' ')[1][1:]
    assert validate_all_conditions(fastqc_installed_version,
                                   version), 'fastqc version is not correct: {} incompatible with {}. Get ' \
                                             'installation instructions here: {}'.format(
        fastqc_installed_version, version_txt, dependencies_installation_page)


def check_cutadapt_exists():
    version, version_txt = dependencies['cutadapt'], dependencies_txt['cutadapt']
    cutadapt_installed_version = run_cmd(['cutadapt', '--version'])
    assert 'cutadapt: not found' not in cutadapt_installed_version, 'cutadapt is not installed'
    assert validate_all_conditions(cutadapt_installed_version,
                                   version), 'cutadapt version is not correct: {} incompatible with {}. Get ' \
                                             'installation instructions here: {}'.format(
        cutadapt_installed_version, version_txt, dependencies_installation_page)


def check_bowtie2_exists():
    version, version_txt = dependencies['bowtie2'], dependencies_txt['bowtie2']
    bowtie2_installed_version = run_cmd(['bowtie2', '--version']).split('\n')[0].split(' ')
    assert bowtie2_installed_version[1] == 'version', 'bowtie2 is not installed'
    assert validate_all_conditions(bowtie2_installed_version[2],
                                   version), 'bowtie2 version is not correct {} incompatible with {}. Get ' \
                                             'installation instructions here: {}'.format(
        bowtie2_installed_version[2], version_txt, dependencies_installation_page)


def check_samtools_exists():
    version, version_txt = dependencies['samtools'], dependencies_txt['samtools']
    samtools_installed_version = run_cmd(['samtools', '--version'])
    assert samtools_installed_version.startswith('samtools'), 'samtools is not installed'
    samtools_installed_version = samtools_installed_version.split('\n')[0].split(' ')[1]
    assert validate_all_conditions(samtools_installed_version,
                                   version), 'samtools version is not correct: {} incompatible with {}. Get ' \
                                             'installation instructions here: {}'.format(
        samtools_installed_version, version_txt, dependencies_installation_page)


def check_rnastructure_exists(path=None):
    version, version_txt = dependencies['rnastructure'], dependencies_txt['rnastructure']
    if path and path is not None:
        cmd = [os.path.join(path, 'Fold'), '--version']
    else:
        cmd = ['Fold', '--version']
    rnastructure_installed_version = run_cmd(cmd)
    assert rnastructure_installed_version.startswith(
        'Fold: Version '), 'RNAstructure is not found. If you are using a Mac and the "Fold" command runs in your ' \
                           'terminal, try giving DREEM the RNAstructure absolute path through the argument ' \
                           '"rnastructure_path" (ex: rnastructure_path="/Users/username/RNAstructure/exe").'
    assert validate_all_conditions(rnastructure_installed_version.split(' ')[2],
                                   version), 'RNAstructure version is not correct: {} incompatible with {}. Get ' \
                                             'installation instructions here: {}'.format(
        rnastructure_installed_version.split(' ')[2], version_txt, dependencies_installation_page)
