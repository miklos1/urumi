from __future__ import absolute_import, print_function, division
from six import iteritems
from six.moves import map, range, zip

from functools import reduce
import itertools
import logging
import signal
from time import time

import numpy

import ufl
import ffc


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
                     TestFunction, det, diff, dot, exp, grad, inner,
                     ln, tr, variable)

    # Define some random parameters
    lamda = Constant(mesh)
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
    def isochoric(F):
        C = F.T*F

        I_1 = tr(C)
        I4_f = dot(e_f, C*e_f)
        I4_s = dot(e_s, C*e_s)
        I8_fs = dot(e_f, C*e_s)

        def cutoff(x):
            return 1.0/(1.0 + exp(-(x - 1.0)*30.0))

        def scaled_exp(a0, a1, argument):
            return a0/(2.0*a1)*(exp(b*argument) - 1)

        E_1 = scaled_exp(a, b, I_1 - 3.)

        E_f = cutoff(I4_f)*scaled_exp(a_f, b_f, (I4_f - 1.)**2)
        E_s = cutoff(I4_s)*scaled_exp(a_s, b_s, (I4_s - 1.)**2)
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
    Fbar = J**(-1.0/3.0)*F

    # Define energy
    E_volumetric = lamda*0.5*ln(J)**2
    psi = isochoric(Fbar) + E_volumetric

    # Find first Piola-Kirchhoff tensor
    P = diff(psi, F)

    # Define the variational formulation
    it = inner(P, grad(v))

    P = ufl.FunctionSpace(mesh, ufl.VectorElement('P', mesh.ufl_cell(), p))
    f = [ufl.Coefficient(P) for _ in range(nf)]
    return ufl.derivative(reduce(ufl.inner,
                                 list(map(ufl.div, f)) + [it])*ufl.dx, u)


forms = [helmholtz, elasticity, hyperelasticity, holzapfel_ogden]


def tsfc_compile_form(form, parameters=None):
    import tsfc
    kernel, = tsfc.compile_form(form, parameters=parameters)
    # FFC generates C strings, TSFC generates a COFFEE AST.
    # Convert COFFEE AST to C string for fairness.
    return kernel.ast.gencode()


def ffc_compile_form(form, parameters=None):
    # Use jit=True to disable element compilation.
    # Makes time comparison fairer.
    if parameters is None:
        parameters = ffc.default_jit_parameters()
    else:
        _ = ffc.default_jit_parameters()
        _.update(parameters)
        parameters = _

    result = ffc.compile_form([form], parameters=parameters, jit=True)
    assert isinstance(result, tuple)
    assert 2 <= len(result) <= 3
    code_h = result[0]
    code_c = result[1]
    # Some versions return a third result which we do not care about.
    return code_c


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

    return ffc.compile_form([form], parameters=parameters)


compilers = {
    'tensor': lambda form: ffc_compile_form(
        form,
        parameters={'representation': 'tensor', 'quadrature_degree': 4}
    ),
    'quadrature': lambda form: ffc_compile_form(
        form,
        parameters={'representation': 'quadrature', 'quadrature_degree': 4}
    ),
    # 'ffc-nonaffine': lambda form: ffc_nonaffine_compile_form(
    #     form,
    #     parameters={'quadrature_degree': 4}
    # ),
    'uflacs': lambda form: ffc_compile_form(
        form,
        parameters={'representation': 'uflacs', 'quadrature_degree': 4}
    ),
    # 'tsfcrepr': lambda form: ffc_compile_form(
    #     form,
    #     parameters={'representation': 'tsfc', 'quadrature_degree': 4}
    # ),
    'tsfc-quick': lambda form: tsfc_compile_form(
        form,
        parameters={'unroll_indexsum': False, 'quadrature_degree': 4}
    ),
    'tsfc-default': lambda form: tsfc_compile_form(
        form,
        parameters={'quadrature_degree': 4}
    ),
}


class TimeoutError(Exception):
    pass


def handler(signum, frame):
    assert signum == signal.SIGALRM
    raise TimeoutError


# Set the signal handler
signal.signal(signal.SIGALRM, handler)


def time_limit(seconds):
    def decorator(func):
        def decorated(*args, **kwargs):
            try:
                signal.alarm(seconds)
                return func(*args, **kwargs)
            finally:
                signal.alarm(0)      # Disable the alarm
        return decorated
    return decorator


def run():
    for domain, form, (name, compiler) in itertools.product(domains, forms, iteritems(compilers)):
        if form in [hyperelasticity, holzapfel_ogden] and name == 'tensor':
            continue
        for nf in range(1):  # 4
            for degree in range(1, 5):
                f = form(domain, degree, degree, nf)

                @time_limit(100)
                def measure():
                    number = 10
                    ends = [None] * number
                    start = time()
                    for i in range(number):
                        compiler(f)
                        ends[i] = time()
                    times = numpy.diff([start] + ends)
                    times.sort()
                    return times[:3].mean()  # average of the three best

                def measure_once():
                    start = time()
                    compiler(f)
                    return time() - start

                try:
                    t = measure()
                except TimeoutError:
                    t = 99.9999
                print(form.__name__, 'dim:', domain.ufl_cell().topological_dimension(), 'degree:', degree, 'nf:', nf, name, t)


if __name__ == "__main__":
    run()
