#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import inspect
import json
import traceback

from sasmodels.core import load_model_info, list_models
from sasmodels.modelinfo import ModelInfo, Parameter
from sasmodels import generate as gen

# Get the current line number so that we can tell the C-compiler
# where to find the broken source code.  Add 2 since that is
# where we use the #line preprocessor command.
LINE = inspect.getlineno(sys._getframe()) + 2
TEMPLATE = """
#line %(LINE)d "webgen.py"
// emscripten uses symbols constant and global which we have redefined
// for OpenCL, so undefine them back to nothing.
#undef constant
#undef global
#include <emscripten/bind.h>
using namespace emscripten;
std::string get_default_params() {
	return "%(INIT_1D)s";
}

typedef std::vector<double> column; // return values holder
std::vector<column>
calculate(%(DECL_1D)s)
{
    std::vector<double> result(q.size());
    for (int i=0; i < q.size(); i++) {
        result[i] = Iq(%(CALL_1D)s);
    }
    return std::vector<column> {result};
}

EMSCRIPTEN_BINDINGS(my_module) {
    register_vector<double>("VectorDouble");
    register_vector<column>("VectorColumn");
    emscripten::function("calculate", &calculate);
    emscripten::function("get_default_params", &get_default_params);
};
"""

def build_cpp(model_info):
    # type: (ModelInfo) -> str
    """
    Like sasmodels.generate, bout outputs an emscripten wrapper instead.
    """
    # Make parameters for q, qx, qy so that we can use them in declarations
    q, qx, qy = [Parameter(name=v) for v in ('q', 'qx', 'qy')]

    partable = model_info.parameters
    pars_1d = [q] + partable.iq_parameters

    # Build parameter and model info
    # p.limits not included because json doesn't handle +/-inf
    init_1d = [[p.id, p.default, p.units, p.description]
               for p in partable.iq_parameters]
    json_info = json.dumps(init_1d).replace('"', r'\"')
    # Note: could construct the html using the following:
    #     make_html(model_info)
    # with the model graph inserted in the usual place.

    # Note: include q in parameter list to simplify code in the case of
    # the porod model which has no parameters.
    decl_1d = ['double %s' % p.id for p in partable.iq_parameters] + [
        "std::vector<double> q",
        "std::vector<double> dq",    # ignored
        "std::vector<double> meanq", # ignored
    ]
    call_1d = ["q[i]"] + [p.id for p in partable.iq_parameters]
    substitutions = {
        'LINE': LINE,
        'INIT_1D': json_info,
        'DECL_1D': ", ".join(decl_1d),
        'CALL_1D': ", ".join(call_1d),
    }
    driver = TEMPLATE % substitutions

    # ... copied from sasmodels.generate.make_source ...
    kernel_header = gen.load_template('kernel_header.c')
    user_code = [(f, open(f).read()) for f in gen.model_sources(model_info)]
    source = []
    gen._add_source(source, *kernel_header)
    for path, code in user_code:
        gen._add_source(source, code, path)

    # Generate form_volume function, etc. from body only
    if isinstance(model_info.form_volume, str):
        pars = partable.form_volume_parameters
        source.append(gen._gen_fn('form_volume', pars, model_info.form_volume,
                                  model_info.filename, model_info._form_volume_line))
    if isinstance(model_info.Iq, str):
        pars = [q] + partable.iq_parameters
        source.append(gen._gen_fn('Iq', pars, model_info.Iq,
                                  model_info.filename, model_info._Iq_line))
    if isinstance(model_info.Iqxy, str):
        pars = [qx, qy] + partable.iqxy_parameters
        source.append(gen._gen_fn('Iqxy', pars, model_info.Iqxy,
                                  model_info.filename, model_info._Iqxy_line))
    # ... done copy ...

    source.append(driver)
    return gen.convert_type("\n".join(source), gen.F64)

PLOT_INSERT = """
.. raw: html

  <div id="top" class="ui-layout-north">Hexagonal</div>
  <div id="center" class="ui-layout-center">
    <div id="plotdiv"></div>
    <div id="plot_controls">
      <label>x-axis
        <select id="xtransform">
          <option>log</option>
          <option>linear</option>
        </select>
      </label>
      <label>y-axis
        <select id="ytransform">
          <option>log</option>
          <option>linear</option>
        </select>
      </label>
      <button id="download_svg">get svg</button>
      <button id="export_model">export</button>
    </div>
  </div>
  <div id="left" class="ui-layout-west">
  </div>

"""
SPLIT_TITLE = '===\n'
def build_html(model_info):
    rst = gen.make_rst(model_info)
    title, body = rst.split(SPLIT_TITLE, 1)


EMCC = "emcc --bind --memory-init-file 0 -s DISABLE_EXCEPTION_CATCHING=0 -O3 -o html/model/%s.js %s"
def compile_model(name):
    # type: (str) -> str
    """
    Return model category, or "unknown" if category is not specified.
    """
    model_info = load_model_info(name)
    model_file = "/tmp/%s.cpp" % model_info.id
    with open(model_file, "w") as fid:
        fid.write(build_cpp(model_info))
    status = os.system(EMCC % (model_info.id, model_file))
    if status == 0:
        return model_info.category.split(':')[-1] if model_info.category else "unknown"
    else:
        return None

MODEL_PATH = "html/model"
CATEGORY_FILE = MODEL_PATH + "/category.json"
def main():
    if len(sys.argv) < 2:
        print("usage: ./webgen.py [model|all]")
        sys.exit(0)

    if not os.path.exists('webgen.py'):
        raise RuntimeError("Must run from the source directory because I'm lazy")

    # Make sure the output directory exists
    if not os.path.exists(MODEL_PATH):
        os.makedirs(MODEL_PATH)

    if os.path.exists(CATEGORY_FILE):
        with open(CATEGORY_FILE) as fid:
            categories = json.load(fid)
    else:
        categories = {}

    model_list = list_models(kind="c") if sys.argv[1] == "all" else sys.argv[1:]
    for name in model_list:
        print("compiling %s..."%name)
        # TODO: maybe continue after error?
        category = compile_model(name)
        if category:
            categories.setdefault(category, []).append(name)
    with open(CATEGORY_FILE, "w") as fid:
        json.dump(categories, fid)

if __name__ == "__main__":
    main()

