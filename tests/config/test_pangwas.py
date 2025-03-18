#!/usr/bin/env python3

import parametrize_from_file as pff
import pangwas
import os
import pytest
import copy
import hashlib
import shutil

from io import StringIO
from contextlib import redirect_stdout


# -----------------------------------------------------------------------------
# Fixtures (from CLI Arguments)

@pytest.fixture(scope="session")
def validate(pytestconfig):
    return pytestconfig.getoption("validate")

@pytest.fixture(scope="session")
def generate(pytestconfig):
    return pytestconfig.getoption("generate")

@pytest.fixture(scope="session")
def safe(pytestconfig):
    return pytestconfig.getoption("safe")

@pytest.fixture(scope="session")
def tags(pytestconfig):
    return pytestconfig.getoption("tags")

@pytest.fixture(scope="session")
def config(pytestconfig):
    return pytestconfig.getoption("config")

# -----------------------------------------------------------------------------
# Tests

@pff
def dev(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def utility(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def extract(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def collect(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def cluster(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def defrag(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def summarize(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def align(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def structural(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def snps(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def presence_absence(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def combine(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def tree(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def root_tree(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def binarize(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def table_to_rtab(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def vcf_to_rtab(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def gwas_lineage(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def gwas_phenotype(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def heatmap(name, params, validate, generate, safe, tags):
    run(**locals())

@pff
def manhattan(name, params, validate, generate, safe, tags):
    run(**locals())

# Annotate is at the end because it's the slowest, and we don't
# use it's output in subsequent steps
@pff
def annotate(name, params, validate, generate, safe, tags):
    run(**locals())

# -----------------------------------------------------------------------------
# Functions

def extract_format_vars(text):
    """
    Extract string format variables from text

    text = '{output_dir}/sample1/{data_dir}'
    extract_format_vars(text)
    ["output_dir", "data_dir"]
    """
    if type(text) != str: return []
    vars = []
    in_var = False
    var = ""
    for c in text:
        if c == "}":
            if in_var == True and var != "":
                vars.append(var)
            in_var = False
            var = ""
        elif c == "{":
            # Ignore nested format strings
            in_var = True if in_var == False else False
        elif in_var:
            var += c
    return vars

def dynamic_format(config:dict, vars:dict, allow_missing:bool=True):
    """
    Apply dynamic string formatting of vars to config.

    :param config: Data dictionary to format
    :param vars:   Variables that hold new string values.
    """
    config = copy.deepcopy(config)
    for k,v in config.items():
        # For dictionary values, start recursion
        if type(v) == dict:        
            v,vars = dynamic_format(config[k], vars, allow_missing)
            config[k] = v
            continue
        # Otherwise, convert all params to a list to iterate through
        # We'll revert back to non-list at the end
        v_is_list = type(v) == list
        values = [v] if type(v) != list else v
        for v_i,v in enumerate(values):
            if type(v) != str: continue
            format_dict = {}
            # Catch if var is already in vars, raise error?
            for var in extract_format_vars(v):
                if var in vars:
                    format_val = vars[var]
                # If it's not a dynamic var, leave it be                    
                elif allow_missing:
                    format_val = "{" + var + "}" 
                else:
                    raise Exception(f"Undefined string formatter {{{var}}} in {k}={v}")
                format_dict[var] = format_val
            if len(format_dict) > 0:
                v = v.format(**format_dict)
            values[v_i] = v
        v = values if v_is_list else values[0]
        vars[k] = v
        config[k] = v

    return (config, vars)

def compare(file_paths:str, validate:bool, generate:bool):
    from distutils.dir_util import copy_tree
    for observed in file_paths:
        expected = observed.replace("observed", "expected")
        if os.path.isdir(observed):
            for file_name in os.listdir(observed):
                file_path = os.path.join(observed, file_name)
                compare([file_path], validate, generate)
        else:
            if generate:
                pangwas.check_output_dir(os.path.dirname(expected))
                if observed != expected:
                    shutil.copyfile(observed, expected)
            if validate:
                print(f"observed: {observed}, expected: {expected}, {sha256sum(observed)[0:8]}, {sha256sum(expected)[0:8]}")                
                assert sha256sum(observed) == sha256sum(expected)


def run_str(s:str, safe:bool=True, vars:dict={}):

    result = None

    # In safe mode, we disable all the builtins
    # like import.
    g = {"__builtins__": {}} if safe else {}
    l = vars

    # Attempt 1: Direct eval
    if result == None:
        #print(f"\tdirect eval: eval({s})")
        try:
            result = eval(f"{s}", g, l)
            #print(f"\tdirect eval result: {result}")
            if type(result) == set:
                result = result.pop()
                #print(f"\tdirect eval set result: {result}")
        except Exception as e:
            #print(f"\tdirect eval failure: {e}")
            pass

    # Attempt 2: exec
    if result == None:
        #print(f"\texec: exec({s})")
        try:
            f = StringIO()
            with redirect_stdout(f):
                exec(f"{s}", g, l)
            result = f.getvalue()
            #print(f"\texec: {result}")
        except Exception as e:
            #print(f"\texec failure: {e}")
            pass
    
    # Fallback: str
    if result == None:
        #print(f"\tstr eval: eval(\"{s}\")")
        result = eval(f"\"{s}\"", g, l)
        #print(f"\tstr eval result: {result}")

    return result

def sha256sum(file_path):
    with open(file_path, 'rb') as infile:
        bytes = infile.read()
        return hashlib.sha256(bytes).hexdigest()

def run(
        name: str,
        params: dict,
        validate:bool=False,
        generate:bool=False,
        safe:bool=True,
        tags:list=[]
    ):

    # Tag filtering
    run_tags = params["tags"] if "tags" in params else []

    if len(tags) > 0:
        exclude = True        
        for tag in tags:
            if tag in run_tags:
                exclude = False
    else:
        exclude = False

    if exclude:
        return

    # Mandatory args
    assert "function" in params
    func =  eval(f"{params['function']}")

    # Check for param keys
    args          = params["args"] if "args" in params else {}
    output        = params["output"] if "output" in params else {}
    variables     = params["variables"] if "variables" in params else None
    error_message = params["error_message"] if "error_message" in params else None
    result        = params["result"] if "result" in params else None

    run_data = [{
        "args": args,
        "output": output,
        "error_message": error_message,
        "function": func,
        "result": result,
    }]

    if type(variables) == dict:
        vars = list(variables.keys())
        for var in vars:
            #print(f"var: {var}")
            values = variables[var]
            # Dict not implemented, skip it
            if type(values) == dict:
                continue
            elif type(values) == str:
                values_run = run_str(values, safe)
                if type(values_run) == list:
                    values = values_run
            values = [values] if type(values) != list else values
            run_data_expanded = []
            for data in run_data:
                for val in values:
                    # In safe mode, we disable all the builtins
                    # like import.
                    g = {"__builtins__": {}} if safe else {}
                    l = locals()
                    result = run_str(s=val, safe=safe, vars={"params": data})
                    exec(f"{var} = '{result}'", g, l)
                    #print(f"\tdata_pre:  {data}")
                    v_data, _vars = dynamic_format(data, vars=l)
                    #print(f"\tdata_post: {v_data}")
                    run_data_expanded.append(v_data)
            run_data = run_data_expanded


    # Final pass, check that no more missing variables are present
    for i,data in enumerate(run_data):
        data, _vars = dynamic_format(data, vars=locals(), allow_missing=False)
        run_data[i] = data

    reproducible_hashes = []
    for data in run_data:
        #print("data:", data)
        args, error_message = data["args"], data["error_message"]
        output = list(data["output"].values())
        func = data["function"]

        # Special handling for annotate, because we don't have a
        # a bakta db available through CI
        if func == pangwas.annotate and "db" in args and not os.path.exists(args["db"]):
            return

        if error_message == None:
            result = func(**args)
            if data["result"] != None:
                assert str(result) == data["result"]
        else:
            try:
                result = func(**args)
                assert False
            except Exception as e:
                assert error_message in f"{e}"

        # Check that all output files exist
        for file_path in output:
            assert os.path.exists(file_path)

        # Hash files for reproducible check
        if "reproducible" in run_tags:
            if len(reproducible_hashes) == 0:
                reproducible_hashes = [
                    {"hashes": set(), "file_paths": []} 
                    for i in range(0,len(output))
                ]
            for i,file_path in enumerate(output):
                hash = sha256sum(file_path)
                reproducible_hashes[i]["hashes"].add(hash)
                reproducible_hashes[i]["file_paths"].append(file_path)
        # We only do observed vs expected comparison on non-reproducible runs
        else:
            compare(output, validate, generate)

    # Final check for reproducibility
    for data in reproducible_hashes:
        assert len(data["hashes"]) == 1 and len(data["file_paths"]) > 1
