from __future__ import absolute_import, print_function, division
from six import iteritems
from six.moves import map, range

from functools import reduce
import itertools
import logging
import re
import signal
import subprocess
import sys
from time import time

import ufl
try:
    import ffc
except ImportError:
    ffc = None


ffc_logger = logging.getLogger("FFC")
ffc_logger.setLevel(logging.WARNING)


try:
    # Check for the presence of new UFL classes
    ufl.Mesh
    ufl.FunctionSpace

except AttributeError:
    # Old UFL detected: set up a compatibility layer.
    def mesh(element):
        return ufl.Domain(ufl.Coefficient(element))

    def function_space(domain, element):
        return element.reconstruct(domain=domain)

    ufl.Mesh = mesh
    ufl.FunctionSpace = function_space
    ufl.Domain.ufl_cell = ufl.Domain.cell


domains = [ufl.Mesh(ufl.VectorElement('P', ufl.triangle, 1)),
           ufl.Mesh(ufl.VectorElement('P', ufl.tetrahedron, 1))]


def helmholtz(domain, q, p, nf=0):
    # Based on https://github.com/firedrakeproject/firedrake-bench/blob/experiments/forms/firedrake_forms.py
    V = ufl.FunctionSpace(domain, ufl.FiniteElement('P', domain.ufl_cell(), q))
    P = ufl.FunctionSpace(domain, ufl.FiniteElement('P', domain.ufl_cell(), p))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    f = [ufl.Coefficient(P) for _ in range(nf)]
    it = ufl.dot(ufl.grad(v), ufl.grad(u)) + 1.0*v*u
    return reduce(ufl.inner, f + [it])*ufl.dx


def elasticity(domain, q, p, nf=0):
    # Based on https://github.com/firedrakeproject/firedrake-bench/blob/experiments/forms/firedrake_forms.py
    V = ufl.FunctionSpace(domain, ufl.VectorElement('P', domain.ufl_cell(), q))
    P = ufl.FunctionSpace(domain, ufl.FiniteElement('P', domain.ufl_cell(), p))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    def eps(v):
        return ufl.grad(v) + ufl.transpose(ufl.grad(v))
    it = 0.25*ufl.inner(eps(v), eps(u))
    f = [ufl.Coefficient(P) for _ in range(nf)]
    return reduce(ufl.inner, f + [it])*ufl.dx


def hyperelasticity(domain, q, p, nf=0):
    # Based on https://github.com/firedrakeproject/firedrake-bench/blob/experiments/forms/firedrake_forms.py
    V = ufl.FunctionSpace(domain, ufl.VectorElement('P', domain.ufl_cell(), q))
    P = ufl.FunctionSpace(domain, ufl.VectorElement('P', domain.ufl_cell(), p))
    v = ufl.TestFunction(V)
    du = ufl.TrialFunction(V)  # Incremental displacement
    u = ufl.Coefficient(V)     # Displacement from previous iteration
    B = ufl.Coefficient(V)     # Body force per unit mass
    # Kinematics
    I = ufl.Identity(domain.ufl_cell().topological_dimension())
    F = I + ufl.grad(u)        # Deformation gradient
    C = F.T*F                  # Right Cauchy-Green tensor
    E = (C - I)/2              # Euler-Lagrange strain tensor
    E = ufl.variable(E)
    # Material constants
    mu = ufl.Constant(domain)  # Lame's constants
    lmbda = ufl.Constant(domain)
    # Strain energy function (material model)
    psi = lmbda/2*(ufl.tr(E)**2) + mu*ufl.tr(E*E)
    S = ufl.diff(psi, E)       # Second Piola-Kirchhoff stress tensor
    PK = F*S                   # First Piola-Kirchoff stress tensor
    # Variational problem
    it = ufl.inner(PK, ufl.grad(v)) - ufl.inner(B, v)
    f = [ufl.Coefficient(P) for _ in range(nf)]
    return ufl.derivative(reduce(ufl.inner,
                                 list(map(ufl.div, f)) + [it])*ufl.dx, u, du)


def holzapfel_ogden(mesh, q, p, nf=0):
    # Based on https://gist.github.com/meg-simula/3ab3fd63264c8cf1912b
    #
    # Original credit note:
    # "Original implementation by Gabriel Balaban,
    # modified by Marie E. Rognes"

    from ufl import (Constant, VectorConstant, Identity, Coefficient,
                     TestFunction, conditional, det, diff, dot, exp,
                     grad, inner, tr, variable)

    # Define some random parameters
    a = Constant(mesh)
    b = Constant(mesh)
    a_s = Constant(mesh)
    b_s = Constant(mesh)
    a_f = Constant(mesh)
    b_f = Constant(mesh)
    a_fs = Constant(mesh)
    b_fs = Constant(mesh)

    # For more fun, make these general vector fields rather than
    # constants:
    e_s = VectorConstant(mesh)
    e_f = VectorConstant(mesh)

    # Define the isochoric energy contribution
    def isochoric(C):
        I_1 = tr(C)
        I4_f = dot(e_f, C*e_f)
        I4_s = dot(e_s, C*e_s)
        I8_fs = dot(e_s, C*e_f)

        def heaviside(x):
            return conditional(x < 1, 0, 1)

        def scaled_exp(a, b, x):
            return a/(2*b)*(exp(b*x) - 1)

        E_1 = scaled_exp(a, b, I_1 - 3)

        E_f = heaviside(I4_f)*scaled_exp(a_f, b_f, (I4_f - 1)**2)
        E_s = heaviside(I4_s)*scaled_exp(a_s, b_s, (I4_s - 1)**2)
        E_3 = scaled_exp(a_fs, b_fs, I8_fs**2)

        E = E_1 + E_f + E_s + E_3
        return E

    # Define mesh and function space
    V = ufl.FunctionSpace(mesh, ufl.VectorElement("CG", mesh.ufl_cell(), q))
    u = Coefficient(V)
    v = TestFunction(V)

    # Misc elasticity related tensors and other quantities
    I = Identity(mesh.ufl_cell().topological_dimension())
    F = grad(u) + I
    F = variable(F)
    J = det(F)
    Cbar = J**(-2/3)*F.T*F

    # Define energy
    Psi = isochoric(Cbar)

    # Find first Piola-Kirchhoff tensor
    P = diff(Psi, F)

    # Define the variational formulation
    it = inner(P, grad(v))

    P = ufl.FunctionSpace(mesh, ufl.VectorElement('P', mesh.ufl_cell(), p))
    f = [ufl.Coefficient(P) for _ in range(nf)]
    return ufl.derivative(reduce(ufl.inner,
                                 list(map(ufl.div, f)) + [it])*ufl.dx, u)


forms = [helmholtz, elasticity, hyperelasticity, holzapfel_ogden]


def tsfc_compile_form(form, parameters=None):
    import tsfc

    tic = time()
    # FFC generates C strings, TSFC generates a COFFEE AST.
    # Convert COFFEE AST to C string for fairness.
    kernel, = tsfc.compile_form(form, parameters=parameters)
    code = kernel.ast.gencode()
    T1 = time() - tic

    print('#include <math.h>\n',
          '#include <string.h>\n',
          code.replace('static inline', '', 1),
          file=open("Form.c", 'w'))

    tic = time()
    subprocess.check_call(["cc", "-pipe", "-c", "-O2", "-std=c99", "Form.c"])
    T2 = time() - tic

    inst = int(subprocess.check_output("objdump -d Form.o | wc -l", shell=True))
    return T1 + T2, inst


def ffc_compile_form(form, parameters=None):
    # Use jit=True to disable element compilation.
    # Makes time comparison fairer.
    if parameters is None:
        parameters = ffc.default_jit_parameters()
    else:
        _ = ffc.default_jit_parameters()
        _.update(parameters)
        parameters = _

    tic = time()
    result = ffc.compile_form([form], parameters=parameters, jit=True)
    T1 = time() - tic

    # Some versions return a third result which we do not care about.
    assert isinstance(result, tuple)
    assert 2 <= len(result) <= 3
    code_h = result[0]  # noqa: F841
    code_c = result[1]
    # lines = ['#include <math.h>', '#include <string.h>']  # TSFC representation
    # lines = ['#include <algorithm>']  # UFLACS representation
    # lines = ['#include <algorithm>', "#include <ufc_geometry.h>"]  # Legacy representations
    lines = ["#include <cstring>", "#include <algorithm>", "#include <ufc_geometry.h>"]
    it = iter(code_c.split('\n'))
    for line in it:
        if '::tabulate_tensor' in line:
            lines.append(re.sub(' .*::', ' ', line))
            for line in it:
                lines.append(re.sub('\) const', ')', line))
                if line == '}':
                    break
            else:
                assert False
    code_c = '\n'.join(lines)
    print(code_c, file=open("Form.cpp", 'w'))

    tic = time()
    subprocess.check_call(["c++", "-pipe", "-c", "-O2", "-std=c++11", "-I/home/mh1714/ffc-tsfc/src/ffc/ffc/backends/ufc/", "Form.cpp"])
    T2 = time() - tic

    inst = int(subprocess.check_output("objdump -d Form.o | wc -l", shell=True))
    return T1 + T2, inst


def ffc_nonaffine_compile_form(form, parameters=None):
    if parameters is None:
        parameters = ffc.default_parameters()
    else:
        _ = ffc.default_parameters()
        _.update(parameters)
        parameters = _

    # Firedrake settings
    parameters["write_file"] = False
    parameters["format"] = 'pyop2'
    parameters["representation"] = 'quadrature'
    parameters["pyop2-ir"] = True

    tic = time()
    kernel, = ffc.compile_form([form], parameters=parameters)
    code = kernel.gencode()
    T1 = time() - tic

    print('#include <math.h>\n',
          code.replace('static inline', '', 1),
          file=open("Form.c", 'w'))

    tic = time()
    subprocess.check_call(["cc", "-pipe", "-c", "-O2", "-std=c99", "Form.c"])
    T2 = time() - tic

    inst = int(subprocess.check_output("objdump -d Form.o | wc -l", shell=True))
    return T1 + T2, inst


compilers_current = {
    'quadrature': lambda form: ffc_compile_form(
        form,
        parameters={'representation': 'quadrature', 'optimize': False, 'quadrature_degree': 4}
    ),
    'uflacs': lambda form: ffc_compile_form(
        form,
        parameters={'representation': 'uflacs', 'optimize': False, 'quadrature_degree': 4}
    ),
    'tsfc-default': lambda form: tsfc_compile_form(
        form,
        parameters={'quadrature_degree': 4}
    ),
    'tsfcrepr': lambda form: ffc_compile_form(
        form,
        parameters={'representation': 'tsfc', 'optimize': False, 'quadrature_degree': 4}
    ),
}


compilers_obsolete = {
    'ffc-nonaffine': lambda form: ffc_nonaffine_compile_form(
        form,
        parameters={'quadrature_degree': 4}
    ),
}


compilers_unused = {
    'tensor': lambda form: ffc_compile_form(
        form,
        parameters={'representation': 'tensor', 'quadrature_degree': 4}
    ),
    'tsfc-quick': lambda form: tsfc_compile_form(
        form,
        parameters={'unroll_indexsum': False, 'quadrature_degree': 4}
    ),
}


class TimeoutError(Exception):
    pass


def handler(signum, frame):
    assert signum == signal.SIGALRM
    raise TimeoutError


# Set the signal handler
signal.signal(signal.SIGALRM, handler)


def time_out(seconds):
    def decorator(func):
        def decorated(*args, **kwargs):
            try:
                signal.alarm(seconds)
                return func(*args, **kwargs)
            finally:
                signal.alarm(0)      # Disable the alarm
        return decorated
    return decorator


def run(compilers, time_limit):
    for domain, form, (name, compiler) in itertools.product(domains, forms, iteritems(compilers)):
        if name == 'tensor' and form in [hyperelasticity, holzapfel_ogden]:
            # Fails due to missing features
            continue

        if name == 'ffc-nonaffine' and form == holzapfel_ogden and domain.ufl_cell().geometric_dimension() > 2:
            # Takes "infinite" time to complete
            continue

        nf = 0
        degree = 2

        f = form(domain, degree, degree, nf)
        number = 18

        @time_out(time_limit)
        def measure():
            ts = [None] * number
            ns = set()
            for i in range(number):
                t, n = compiler(f)
                ts[i] = t
                ns.add(n)
            n, = ns
            return ts, n

        def measure_once():
            t, n = compiler(f)
            return [t], n

        try:
            times, size = measure()

            for t in times:
                print(form.__name__, domain.ufl_cell().topological_dimension(), name, t, size)
        except TimeoutError:
            print('Measurement timed out:', form.__name__, domain.ufl_cell().topological_dimension(), name,
                  file=sys.stderr)


if __name__ == "__main__":
    _, mode = sys.argv
    if mode == 'current':
        run(compilers_current, 8640)
    elif mode == 'ffc-bendy':
        run(compilers_obsolete, 18000)
    else:
        assert False, "unrecognised mode"
