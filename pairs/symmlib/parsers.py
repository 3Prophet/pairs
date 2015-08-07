import os
import re

def line_yielder(fpath):
    with open(fpath,'r') as fobj:
        for line in fobj:
            yield line

def line_splitter(fpath):
    for line in line_yielder(fpath):
        splitted_line =  line.strip().split()
        yield splitted_line

def filtered_columns(fpath):
    path

def assigner(keys):
    pass

def gen_open(fpaths):
    for fpath in fpaths: 
        yield open(fpath, 'r')

def gen_lines(fobjs):
    for fobj in fobjs:
        for line in fobj:
            yield line

def gen_split(lines):
    for line in lines:
        yield line.strip().split()

def gen_zip(lines, keys):
    for line in lines:
        yield dict(zip(keys,line))

def gen_selector(dicts, keys):
    for d in dicts:
        selected = {}
        for k in keys:
            selected[k] = d[k]
        yield selected

def split_parse(fpath, keys_inp, keys_out):
    fobjs = gen_open((fpath,))
    lines = gen_lines(fobjs)
    splited_lines = gen_split(lines)
    zipped_lines = gen_zip(splited_lines,keys_inp)
    selected_lines = gen_selector(zipped_lines, keys_out)
    for line in selected_lines:
        yield line

if __name__ == "__main__":
    dpath = "/Users/dima/DiffuseScattering/SNBL/beta_bragg_XDS/shelx/"
    fname1 = "beta.fcf"
    keys_inp = ('h','k','l','F_o','F_c','phi')
    keys_out = ('h','k','l','F_o')
    p = multiprocessing.Process(target = split_parse,
                    args = (os.path.join(dpath, fname1), keys_inp, keys_out),
                    )
    p.daemon =True
    p.start()

