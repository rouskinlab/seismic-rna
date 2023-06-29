import numpy as np
import datetime
from os.path import exists
import os
import matplotlib.pyplot as plt
from typing import Tuple, List
import shutil
from plotly.validators.scatter.marker import SymbolValidator
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM
from scipy.optimize import curve_fit
import scipy
import inspect

vals = SymbolValidator().values

def Setshape(x):
        vals = SymbolValidator().values
        return vals[3*x]


def make_path(path:str)->str:
    """Create directories until path exists on your computer. Turns the keyword 'date' into today's date.

    Args:
        path: series of directories that you want to create.
    
    Returns:
        Updated path with today's date instead of the keyword 'date'  
    """

    path = os.path.normpath(path)
    path=path.split(os.sep)
    try:
        path[path.index('date')] = str(datetime.datetime.now())[:10]
    except:
        'No date in path'
    full_path = ''
    for repo in path:
        full_path = full_path + f"{repo}/"
        if not exists(full_path):
            os.mkdir(full_path)
    return full_path



def gini(x:np.array)->float:
    """Returns Gini index

    Args:
        x (np.array): the array you want the Gini index from

    Returns:
        float: Gini index of the input array
    """
    # (Warning: This is a concise implementation, but it is O(n**2)
    # in time and memory, where n = len(x).  *Don't* pass in huge
    # samples!)

    # Mean absolute difference
    mad = np.abs(np.subtract.outer(x, x)).mean()
    # Relative mean absolute difference
    rmad = mad/np.mean(x)
    # Gini coefficient
    g = 0.5 * rmad
    return g

def savefig(file:str, close=True)->None:
    """Save a matplotlib figure and create the directory if it doesn't exists.

    Args:
        file: path+title.
        facecolor: color of the background 
    """

    path = make_path('/'.join(file.split('/')[:-1]))
    plt.savefig(path+file.split('/')[-1], bbox_inches='tight')
    if close:
        # Clear the current axes.
        plt.cla() 
        # Clear the current figure.
        plt.clf() 
        # Closes all the figure windows.
        plt.close('all')   


def define_figure(title:str, xlabel:str, ylabel:str, figsize:Tuple[float, float])->plt.figure:
    """Define title, labels and size of your figure.

    Args:
        title: matplotlib title
        xlabel: matplotlib xlabel
        ylabel: matplotlib ylabel
        figsize: matplotlib figsize
    """

    fig = plt.figure(figsize=figsize)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    return fig

import yaml, os, subprocess
import pandas as pd


class Fit(object):
    def __init__(self) -> None:
        self.legend = ''
    
    def get_legend(self):
        return self.legend

    def predict(self, x, y, model, prefix='', suffix=''):
        fit = self.fit(x,y,model)
        m = eval(model)
        try:
            linreg  = scipy.stats.linregress(y,m(x,*fit))
            self.rvalue = round(linreg.rvalue,5)
        except:
            self.rvalue = 'error'

        self._generate_legend(fit, model, prefix, suffix)
        return np.sort(x), m(np.sort(x),*fit)

    def fit(self, x,y, model):
        fit = curve_fit(eval(model), x, y)[0]
        return fit

    def _generate_legend(self, fit, m, prefix, suffix):
        slice_m = lambda start, stop: ','.join(str(m).split(',')[start:stop])
        first_slice = slice_m(0,len(fit))+','+slice_m(len(fit), len(fit)+1).split(':')[0]
        second_slice = ','.join(m.split(',')[len(fit):])[2:]
        fun_args = [a.strip() for a in str(m).split(',')[1:len(fit)+1]]
        fun_args[-1] = fun_args[-1][0]
        for a,b in zip(fun_args,fit):
            second_slice = second_slice.replace(a.strip(),str(round(b,5)))
        self.legend = prefix+ second_slice + suffix +f'\n R2={self.rvalue}'

                
class RNAstructure(object): 
    """A class to run RNAstructure commands

    Args:
        object (_type_): _description_

    Examples:
        rna = RNAstructure('/Users/ymdt/dreem/RNAstructure/exe')
        rna.fit(sequence='ACCTTTCAGAGCTACGATCGACTAGCTAGCATCGATACAGCGACACAAGCATTTGTAGCATTAGGTCA')
        print("DeltaG + structure:", rna.predict_reference_deltaG())
        print("Ensemble energy:", rna.predict_ensemble_energy())
        print("Partition function:", rna.predict_partition()) 
    """
    def __init__(self, rnastructure_path) -> None:
        self.rnastructure_path = rnastructure_path if rnastructure_path[-1] == '/' else rnastructure_path+'/'

    def fit(self, sequence, reference='reference'):
        self.sequence = sequence
        temp_folder = 'temp/rnastructure/'
        self.__make_temp_folder(temp_folder)
        self.__make_files(temp_folder+'temp')
        self.__create_fasta_file(reference, sequence)

    def predict_ensemble_energy(self):
        cmd = f"{self.rnastructure_path}EnsembleEnergy {self.fasta_file} --DNA --sequence"
        splitted_output = self.__run_command(cmd)[0].split(' ')
        return float(splitted_output[splitted_output.index(f"kcal/mol\n\nEnsemble")-1])

    def predict_partition(self, temperature_k =None):
        cmd = f"{self.rnastructure_path}partition {self.fasta_file} {self.pfs_file} --DNA"
        if temperature_k != None:
            cmd += ' --temperature '+str(temperature_k)
        self.__run_command(cmd)
        self.__run_command(self.rnastructure_path+'ProbabilityPlot '+ self.pfs_file + ' -t '+self.prob_file)
        with open(self.prob_file,"r") as f:
            lines=f.readlines()
            out={'i':[],'j':[],'p':[]}
            for x in range(len(lines)):
                if x>1:
                    ls=lines[x].split("\t")
                    out["i"]+=[int(ls[0])]
                    out["j"]+=[int(ls[1])]
                    out["p"]+=[float(ls[2])]
        return self.__cast_pairing_prob(out)

    def draw(self, savefig, sub_rate = [], dpi=72):
        self.predict_reference_deltaG()
        if len(sub_rate):
            sub_rate = np.array(sub_rate).reshape(-1,1)
            with open(self.s_file, 'w') as f:
                for i in range(len(sub_rate)):
                    f.write(f'{i+1} {sub_rate[i][0] if sub_rate[i][0] != np.nan else -999 } \n')
                f.close()
            
        cmd = self.rnastructure_path+'draw '+self.ct_file + ' --svg '+self.svg_file+' -n -1'
        self.__run_command(cmd if len(sub_rate) == 0 else cmd+' -s '+self.s_file)
        svg_code = self.__read_svg()
        drawing = svg2rlg(self.svg_file)
        renderPM.drawToFile(drawing, self.png_file, fmt="PNG", dpi=dpi)
        shutil.copy(self.png_file, savefig)

    def __read_svg(self):
        with open(self.svg_file, "r") as f:
            svg_code = f.read()
        return svg_code

    def predict_reference_deltaG(self, temperature_k=None):
        # Run RNAstructure
        suffix = ''
        cmd = f"{self.rnastructure_path}Fold {self.fasta_file} {self.ct_file} -d" + suffix
        self.__run_command(cmd)
        assert os.path.getsize(self.ct_file) != 0, f"{self.ct_file} is empty, check that RNAstructure works"
        self.deltaG, self.structure =  self.__extract_deltaG_struct()
        return self.deltaG, self.structure

    def __make_temp_folder(self, temp_folder):
        isExist = os.path.exists(temp_folder)
        if not isExist:
            os.makedirs(temp_folder)
        return temp_folder

    def __make_files(self, temp_prefix='temp/temp'):
        self.pfs_file = f"{temp_prefix}.pfs"
        self.ct_file = f"{temp_prefix}.ct"
        self.dot_file = f"{temp_prefix}_dot.txt"
        self.fasta_file = temp_prefix+'fasta'
        self.prob_file = temp_prefix+'_prob.txt'
        self.svg_file = temp_prefix+'.svg'
        self.png_file = temp_prefix+'.png'
        self.s_file = temp_prefix+'_s.txt'

    def __create_fasta_file(self, reference, sequence):
        # push the ref into a temp file
        temp_fasta = open(self.fasta_file, 'w')
        temp_fasta.write('>'+reference+'\n'+sequence)
        temp_fasta.close()

    # cast the temp file into a dot_bracket structure and extract the attributes
    def __extract_deltaG_struct(self):
        self.__run_command(f"ct2dot {self.ct_file} 1 {self.dot_file}")
        temp_dot = open(self.dot_file, 'r')
        first_line = temp_dot.readline().split()
        # If only dots in the structure, no deltaG 
        out = {}
        if len(first_line) == 4:
            _, _, deltaG, _ = first_line
            deltaG = float(deltaG)
        if len(first_line) == 1:
            deltaG, _ = 'void', first_line[0][1:]

        sequence = temp_dot.readline()[:-1] #  Remove the \n
        structure = temp_dot.readline()[:-1] # Remove the \n
        return deltaG, structure

    def __cast_pairing_prob(self, prob:dict)->list:
        """Computes the pairing probability list.
        Args:
            prob (dict): Ouput of RNAstructure.
            index (int): index of the row that you want the pairing probabilities of.
        Returns:
            list: pairing probability, under the form of a list of probabilities
        """
        # Create local dataframe to play here
        df_loc = pd.DataFrame(prob)
        df_loc = pd.concat((df_loc, df_loc.rename(columns={'i':'j','j':'i'})))
        # Group the probabilities by i and get an ordered list of the pairing probability
        g = df_loc.groupby('i')['p'].agg(lambda row: [pow(10,-float(r)) for r in row])
        g.index = g.index.astype(int)
        g = g.sort_index()
        g['sum_log_p'] = g.apply(lambda row: sum(row))
        return list(g['sum_log_p'])

    def __run_command(self, cmd):
        output, error_msg = None, None
        try:
            output = subprocess.check_output(
                    cmd, shell=True, stderr=subprocess.STDOUT
            ).decode("utf8")
        except subprocess.CalledProcessError as exc:
            error_msg = exc.output.decode("utf8")
        return output, error_msg

def save_plot(func,  *args, **kwargs):
    """Default arguments for saving a plotly figure.

    Args:
        to_html (str, optional): File name to save the figure as a html.
        to_png (str, optional): File name to save the figure as a png.
    """

    def wrapper(*args, **kwargs):
        out = func(*args, **kwargs)
        to_html, to_png = kwargs.pop('to_html', None), kwargs.pop('to_png', None)
        if to_html:
            out['fig'].write_html(to_html)
        if to_png:
            out['fig'].write_image(to_png)
        return out
    wrapper.__name__=func.__name__
    wrapper.__doc__=func.__doc__
    return wrapper   

def extract_args(func):
    """Extract the arguments of a function.

    Args:
        func (function): Function to extract the arguments from.

    Returns:
        list: List of the arguments of the function.
    """
    return inspect.getfullargspec(func).args


def assert_only_one_row(df):
    """Assert that the dataframe has only one row.

    Args:
        df (pd.DataFrame): Dataframe to check.

    Raises:
        ValueError: If the dataframe has more than one row.
    """
    if df.shape[0] > 1:
        raise ValueError("The dataframe has more than one row.")
    
    if df.shape[0] == 0:
        raise ValueError("The dataframe is empty.")
