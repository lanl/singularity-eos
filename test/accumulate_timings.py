#!/usr/bin/env python

import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser(description="Combine outputs from compare_to_eospac")
parser.add_argument('--names',type=str,metavar='N',nargs='+',
                    help='names of materials')
parser.add_argument('--matids',type=int,metavar='M',nargs='+',
                    help='matides of materials')
parser.add_argument('--files',type=str,metavar='F',nargs='+',
                    help='filenames')
parser.add_argument('--sizes',type=int,metavar='S',nargs='+',
                    help='sizes of tests')

class Timings:
    def __init__(self,diff=None,pac=None,cpu=None,gpu=None):
        self.diff = diff
        self.pac = pac
        self.cpu = cpu
        self.gpu = gpu
        self._keys = ['diff','pac','cpu','gpu']
    def keys(self):
        return self._keys
    def __getitem__(self,key):
        if key == 'diff':
            return self.diff
        elif key == 'pac':
            return self.pac
        elif key == 'cpu':
            return self.cpu
        elif key == 'gpu':
            return self.gpu
        else:
            raise ValeError("Key must be one of {}".format(self._keys))

def get_val(line):
    return float(line.split(' = ')[1].lstrip().rstrip())

def parse_file(filename,matids):
    interp = {matid : Timings() for matid in matids}
    root = {matid : Timings() for matid in matids}
    current_matid = None
    root_finding = False
    with open(filename,'r') as f:
        for line in f:
            for matid in matids:
                if '...{}'.format(matid) in line:
                    current_matid = matid
                    root_finding = False
            if current_matid is not None:
                if 'Root finding:' in line:
                    root_finding = True
                if '...L2 difference = ' in line:
                    interp[current_matid].diff = get_val(line)
                if '...EOSPAC time/point (microseconds) = ' in line:
                    if root_finding:
                        root[current_matid].pac = get_val(line)
                    else:
                        interp[current_matid].pac = get_val(line)
                if '...spiner host time/point (microseconds) = ' in line:
                    interp[current_matid].cpu = get_val(line)
                if '...spiner device time/point (microseconds) = ' in line:
                    interp[current_matid].gpu = get_val(line)
                if '...L2 difference for spiner(rho,T) host = ' in line:
                    root[current_matid].diff = get_val(line)
                if '...spiner(rho,T) time/point (microseconds) host = ' in line:
                    root[current_matid].cpu = get_val(line)
                if '...spiner(rho,T) time/point (microseconds) device = ' in line:
                    root[current_matid].gpu = get_val(line)
    return interp,root

if __name__ == "__main__":
    args = parser.parse_args()
    assert (args.matids),"there must be materials"
    assert (args.files),"there must be files"
    assert (len(args.sizes) == len(args.files)),"files must be tied to sizes"
    assert (len(args.matids) == len(args.names)),"names must be tied to matids"

    names = {}
    for name,matid in zip(args.names,args.matids):
        names[matid] = name
    data = np.empty((len(args.sizes),len(args.matids)*2*4+1))
    for s,(size,fnam) in enumerate(zip(args.sizes,args.files)):
        data[s,0] = size
        interp,root = parse_file(fnam,args.matids)
        for i,matid in enumerate(args.matids):
            stride = i*4*2
            for j,key in enumerate(interp[matid].keys()):
                data[s,stride+j+1] = interp[matid][key]
            for j,key in enumerate(root[matid].keys()):
                data[s,stride+1+4+j] = root[matid][key]
    header_md = '| size |'
    header_txt = 'size'
    for matid in args.matids:
        for tpe in ['interp','root']:
            for key in interp[matid].keys():
                header_md += '{}_{}_{} |'.format(matid,tpe,key)
                header_txt += '\t{}_{}_{}'.format(matid,tpe,key)
    np.savetxt('timings.md',data,
               delimiter = ' | ',
               header = header_md,
               comments = '')
    np.savetxt('timings.dat',data,
               header = header_txt,
               comments = '# ')

