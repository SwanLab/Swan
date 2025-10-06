import numpy as np
import numpy.testing as npt
import sys

from nullspace_optimizer.examples.basic_examples import ex00_basic_problem, \
    ex01_basic_problem_2,\
    ex02_basic_problem_3, \
    ex03_basic_problem_degenerate, \
    ex05_multiple_shootings, \
    ex06_multiple_shootings2, \
    ex07_unconstrained, \
    ex08_degenerate

from nullspace_optimizer import nlspace_solve
from nullspace_optimizer.optimizers.nullspace.utils import compute_norm


def test_ex00_basic_problem_1():
    results = ex00_basic_problem.solve_basic_problem()
    actual = results['x'][-1]
    expected = np.asarray([1.82569891, 0.54773544])
    npt.assert_allclose(actual, expected, atol=1.0e-5)

test_ex00_basic_problem_1(), \

def test_ex00_basic_problem_equalized():
    results = ex00_basic_problem.solve_basic_problem_equalized()
    actual = results['x'][-1][:2]
    expected = np.asarray([1.82572629, 0.54772723])
    npt.assert_allclose(actual, expected, atol=3.0e-5, rtol=0)


def test_ex00_basic_problem_restart():
    results = ex00_basic_problem.solve_basic_problem()
    resultsRestart = ex00_basic_problem.restart_basic_problem(results)
    for key in results:
        npt.assert_equal(len(resultsRestart[key]), len(results[key]))
    npt.assert_allclose(resultsRestart['x'][-1], results['x'][-1])


def test_ex01_basic_problem_2():
    results = ex01_basic_problem_2.solve_basic_problem_2()
    actual = results['x'][-1]
    expected = np.asarray([1.5, 1.5])
    npt.assert_allclose(actual, expected, rtol=1.0e-5)


def test_ex01_basic_problem_2_equalized():
    results = ex01_basic_problem_2.solve_basic_problem_2_equalized()
    actual = results['x'][-1][:2]
    expected = np.asarray([1.5, 1.5])
    npt.assert_allclose(actual, expected, rtol=1.0e-5)


def test_ex02_basic_problem_3():
    results = ex02_basic_problem_3.solve_problem_parab()
    actual = results['x'][-1]
    expected = np.asarray([-2.5, 0.5])
    npt.assert_allclose(actual, expected, atol=1.0e-4)


def test_ex02_basic_problem_3_equalized():
    results = ex02_basic_problem_3.solve_problem_parab_equalized()
    actual = results['x'][-1][:2]
    expected = np.asarray([-2.5, 0.5])
    npt.assert_allclose(actual, expected, atol=3.0e-5)


def test_ex02_basic_problem_3_osqp():
    results = ex02_basic_problem_3.solve_problem_parab_osqp()
    actual = results['x'][-1]
    expected = np.asarray([-2.5, 0.5])
    npt.assert_allclose(actual, expected, atol=1.0e-4)

def test_ex03_basic_problem_degenerate():
    results = ex03_basic_problem_degenerate.run_basic_problem()
    actual = results['x'][-1]
    expected = np.asarray([1.82572629, 0.54772723])
    npt.assert_allclose(actual, expected, atol=1.0e-4)


def test_ex03_basic_problem_degenerate_equalized():
    results = ex03_basic_problem_degenerate.run_basic_problem_equalized()
    actual = results['x'][-1][:2]
    expected = np.asarray([1.82572629, 0.54772723])
    npt.assert_allclose(actual, expected, rtol=3.0e-5)


def test_ex05_multiple_shootings():
    results = ex05_multiple_shootings.run_problems()
    expected = np.asarray([-np.sqrt(1 - 0.7**2), -0.7])
    for r in results:
        actual = r['x'][-1]
        npt.assert_allclose(actual, expected, rtol=1e-6)

def test_ex05_multiple_shootings_cvxopt():
    results = ex05_multiple_shootings.run_problems(qp_solver="cvxopt")
    expected = np.asarray([-np.sqrt(1 - 0.7**2), -0.7])
    for r in results:
        actual = r['x'][-1]
        npt.assert_allclose(actual, expected, rtol=1e-6)

def test_ex06_multiple_shootings_2():
    results = ex06_multiple_shootings2.run_problems(dt=0.1)
    solutions = ([0.5, 0.0], [-1.0, -1.2], [-1.0, 0.5],
                 [-1.0, -0.5], [-1.0, 1.0])
    for r, e in zip(results, solutions):
        actual = r['x'][-1]
        expected = np.asarray(e)
        assert np.isclose(actual[0],-1,rtol=1e-6) or np.allclose(actual,[0.5,0.0])


def test_ex07_unconstrained():
    results = ex07_unconstrained.run_problems()
    actual = results['x'][-1]
    expected = np.zeros((2,))
    npt.assert_allclose(actual, expected, atol=5.0e-4)
        

def test_ex08_degenerate():
    results = ex08_degenerate.run_problem()
    actual = results['x'][-1]
    expected = [0,1]
    npt.assert_allclose(actual, expected, atol=5.0e-4)
        
